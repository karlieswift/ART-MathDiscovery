#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, os, re, shutil, subprocess, sys
from pathlib import Path
from common import *
from sixframe_orfs import sixframe

RADIUS=10000; EMAX=1e-4; MINAA=45

def run(cmd,**kw):
    q=subprocess.run(cmd,text=True,capture_output=True,**kw)
    if q.returncode!=0:
        raise RuntimeError('command failed: '+repr(cmd)+'\nSTDOUT\n'+q.stdout[-4000:]+'\nSTDERR\n'+q.stderr[-4000:])
    return q

def load_manifest(here):
    rows=read_tsv(here/'frozen/BASELINE_SEQUENCE_MANIFEST.tsv')
    return {r['case_id']:r for r in rows}

def parse_prodigal_faa_coords(faa:Path):
    # Prodigal protein headers are: >query # start(1-based) # end(1-based) # strand(1/-1) # ...
    coords={}
    for ln in faa.read_text(errors='replace').splitlines():
        if not ln.startswith('>'): continue
        parts=[x.strip() for x in ln[1:].split('#')]
        if len(parts)<4: continue
        q=parts[0].split()[0]
        try:
            start1=int(parts[1]); end1=int(parts[2]); strand='+' if int(parts[3])>=0 else '-'
        except Exception:
            continue
        coords[q]=(start1-1,end1,strand)
    return coords

def parse_tblout(tbl:Path):
    hits=[]
    for ln in tbl.read_text(errors='replace').splitlines():
        if not ln or ln.startswith('#'): continue
        f=ln.split()
        if len(f)<8: continue
        # hmmscan tblout: target=model, query=protein; full seq E-value at col 5 (0-based 4), score col 6 (5)
        model=f[0]; query=f[2]
        try: e=float(f[4]); score=float(f[5])
        except: continue
        if e<=EMAX: hits.append((query,model,e,score))
    # best per protein
    best={}
    for h in hits:
        q=h[0]
        if q not in best or (h[2],-h[3])<(best[q][2],-best[q][3]): best[q]=h
    return list(best.values())

def find_recurrence_proposals(board_inventory_row):
    paths=[Path(x) for x in board_inventory_row.get('proposal_files','').split('|') if x]
    rows=[]
    for p in paths:
        try:
            with p.open(newline='',errors='replace') as f:
                rd=csv.DictReader(f,delimiter='\t')
                for r in rd:
                    txt=' '.join((str(k)+'='+str(v)) for k,v in r.items()).lower()
                    if 'recurr' not in txt: continue
                    # flexible coordinates
                    s=None;e=None;q=None
                    for k,v in r.items():
                        kl=k.lower(); vv=str(v).strip()
                        if s is None and kl in {'start','start0','window_start','window_start0','genomic_start','begin'}:
                            try:s=int(float(vv))
                            except: pass
                        if e is None and kl in {'end','end0','window_end','window_end0','genomic_end','stop'}:
                            try:e=int(float(vv))
                            except: pass
                        if q is None and ('q_recurr' in kl or kl in {'q','qvalue','q_value'}):
                            try:q=float(vv)
                            except: pass
                    if s is not None and e is not None: rows.append((s,e,q if q is not None else 1.0,str(p)))
        except Exception: continue
    return rows

def local_rec(rt0,rt1,rec):
    a=max(0,rt0-RADIUS); b=rt1+RADIUS
    hits=[x for x in rec if x[0]<b and x[1]>a]
    if not hits: return (0,1.0,'')
    best=min(hits,key=lambda x:x[2])
    return (len(hits),best[2],best[3])

def score_case(case, seqrow, boardrow, rvt:Path, outroot:Path, cpu:int):
    cid=case['case_id']; cdir=outroot/cid; cdir.mkdir(parents=True,exist_ok=True)
    fa=Path(seqrow['chosen_fasta'])
    rec=find_recurrence_proposals(boardrow)
    allrows=[]
    # Prodigal
    pdir=cdir/'prodigal'; pdir.mkdir(exist_ok=True)
    gff=pdir/'genes.gff'; faa=pdir/'proteins.faa'; tbl=pdir/'rvt.tblout'; dom=pdir/'rvt.domtblout'
    run(['prodigal','-i',str(fa),'-o',str(gff),'-a',str(faa),'-f','gff','-p','meta','-q'])
    run(['hmmscan','--cpu',str(cpu),'-E',str(EMAX),'--domE',str(EMAX),'--tblout',str(tbl),'--domtblout',str(dom),str(rvt),str(faa)])
    coords=parse_prodigal_faa_coords(faa)
    for q,model,e,score in parse_tblout(tbl):
        if q not in coords: continue
        s,t,strand=coords[q]; nh,rq,rp=local_rec(s,t,rec)
        allrows.append({'case_id':cid,'method':'PRODIGAL_RVT','candidate_id':q,'start0':s,'end0':t,'strand':strand,'hmm_model':model,'hmm_evalue':e,'hmm_score':score,'local_recurrence_hit_n':nh,'best_local_recurrence_q':rq,'recurrence_source':rp})
    # Six-frame
    recs=list(fasta_records(fa)); assert len(recs)==1
    seq=recs[0][1]
    sfdir=cdir/'sixframe'; sfdir.mkdir(exist_ok=True); sffaa=sfdir/'sixframe.faa'; sftbl=sfdir/'rvt.tblout'; sfdom=sfdir/'rvt.domtblout'
    coord={}
    with sffaa.open('w') as f:
        for i,(strand,frame,a,b,s,t,pep) in enumerate(sixframe(seq,MINAA),1):
            q=f'SF{i:08d}'
            coord[q]=(s,t,strand)
            f.write(f'>{q} start0={s} end0={t} strand={strand} frame={frame}\n{pep}\n')
    run(['hmmscan','--cpu',str(cpu),'-E',str(EMAX),'--domE',str(EMAX),'--tblout',str(sftbl),'--domtblout',str(sfdom),str(rvt),str(sffaa)])
    for q,model,e,score in parse_tblout(sftbl):
        if q not in coord: continue
        s,t,strand=coord[q]; nh,rq,rp=local_rec(s,t,rec)
        allrows.append({'case_id':cid,'method':'SIXFRAME_RVT','candidate_id':q,'start0':s,'end0':t,'strand':strand,'hmm_model':model,'hmm_evalue':e,'hmm_score':score,'local_recurrence_hit_n':nh,'best_local_recurrence_q':rq,'recurrence_source':rp})
    fields=['case_id','method','candidate_id','start0','end0','strand','hmm_model','hmm_evalue','hmm_score','local_recurrence_hit_n','best_local_recurrence_q','recurrence_source']
    write_tsv(cdir/'BASELINE_RT_CANDIDATES.tsv',allrows,fields)
    return {'case_id':cid,'n_candidates':len(allrows),'prodigal_n':sum(r['method']=='PRODIGAL_RVT' for r in allrows),'sixframe_n':sum(r['method']=='SIXFRAME_RVT' for r in allrows),'recurrence_proposal_n':len(rec)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--cpu',type=int,default=16); a=ap.parse_args()
    here=Path(__file__).resolve().parents[1]
    gate=json.loads((here/'frozen/BASELINE_B0_GATE.json').read_text())
    if not gate['gate_pass']: raise SystemExit('B0 gate is HOLD; do not score.')
    cases=read_tsv(here/'config/CASES.tsv'); seqmap=load_manifest(here)
    boards={r['case_id']:r for r in read_tsv(here/'frozen/BASELINE_BOARD_INVENTORY.tsv')}
    tool=json.loads((here/'frozen/BASELINE_TOOLCHAIN.json').read_text()); rvt=Path(tool['rvt_path'])
    out=here/'outputs/score_freeze'; out.mkdir(parents=True,exist_ok=True)
    summary=[]
    # CRITICAL: this runner never reads BASELINE_TRUTH_SEAL.tsv
    for c in cases:
        summary.append(score_case(c,seqmap[c['case_id']],boards[c['case_id']],rvt,out,a.cpu))
        print('[SCORE]',summary[-1],flush=True)
    write_tsv(out/'BASELINE_SCORE_SUMMARY.tsv',summary,list(summary[0]))
    hashes=[]
    for p in sorted(out.rglob('*')):
        if p.is_file(): hashes.append({'path':str(p.relative_to(here)),'sha256':sha256_file(p)})
    write_tsv(out/'BASELINE_SCORE_FREEZE_SHA256.tsv',hashes,['path','sha256'])
    gate={'stage':'Stage05F4A-B1-score-freeze','gate':'PASS_STAGE05F4A_B1_PROTEIN_FIRST_BASELINE_SCORES_FROZEN_NO_TRUTH','gate_pass':True,'truth_opened_for_scoring':False,'retuning_authorized':False,'rt_neighborhood_radius_bp':RADIUS,'rvt_evalue':EMAX,'sixframe_min_aa':MINAA}
    safe_json(out/'BASELINE_B1_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
