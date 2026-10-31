#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, subprocess
from pathlib import Path
from common import *
from sixframe_orfs import sixframe

RADIUS=10000; EMAX=1e-4; DOMMAX=1e-4; MINAA=45

def run(cmd):
    q=subprocess.run(cmd,text=True,capture_output=True)
    if q.returncode!=0: raise RuntimeError('command failed: '+repr(cmd)+'\nSTDOUT\n'+q.stdout[-4000:]+'\nSTDERR\n'+q.stderr[-4000:])
    return q

def parse_prodigal_faa_coords(faa:Path):
    coords={}
    for ln in faa.read_text(errors='replace').splitlines():
        if not ln.startswith('>'):continue
        parts=[x.strip() for x in ln[1:].split('#')]
        if len(parts)<4:continue
        q=parts[0].split()[0]
        try:s1=int(parts[1]); e1=int(parts[2]); strand='+' if int(parts[3])>=0 else '-'
        except Exception:continue
        coords[q]=(s1-1,e1,strand)
    return coords

def parse_domtblout(p:Path):
    best={}
    for ln in p.read_text(errors='replace').splitlines():
        if not ln or ln.startswith('#'):continue
        f=ln.split()
        if len(f)<22:continue
        model=f[0]; query=f[3]
        try:full_e=float(f[6]); full_score=float(f[7]); i_e=float(f[12]); dom_score=float(f[13])
        except Exception:continue
        if full_e>EMAX or i_e>DOMMAX:continue
        h=(query,model,full_e,i_e,full_score,dom_score)
        if query not in best or (full_e,i_e,-full_score,-dom_score)<(best[query][2],best[query][3],-best[query][4],-best[query][5]):best[query]=h
    return list(best.values())

def load_rec(here,cid):
    p=here/'frozen/canonical'/cid/'CANONICAL_PROPOSALS.tsv'
    rows=[]
    for r in read_tsv(p):
        if 'RECURRENCE' not in r['family']:continue
        rows.append((int(r['start0']),int(r['end0']),float(r['qvalue']),r['source_path']))
    return rows

def local_rec(s,e,rec):
    a=max(0,s-RADIUS); b=e+RADIUS; hits=[x for x in rec if x[0]<b and x[1]>a]
    if not hits:return 0,1.0,''
    z=min(hits,key=lambda x:x[2]); return len(hits),z[2],z[3]

def union_len(iv):
    if not iv:return 0
    z=sorted(iv); total=0; a,b=z[0]
    for c,d in z[1:]:
        if c<=b:b=max(b,d)
        else:total+=b-a; a,b=c,d
    return total+b-a

def score_case(case,seqrow,rvt,prodigal,hmmscan,outroot,cpu,here):
    cid=case['case_id']; cdir=outroot/cid; cdir.mkdir(parents=True,exist_ok=True)
    fa=Path(seqrow['chosen_fasta']); n=int(seqrow['length_bp']); rec=load_rec(here,cid); rows=[]
    pdir=cdir/'prodigal'; pdir.mkdir(exist_ok=True); gff=pdir/'genes.gff'; faa=pdir/'proteins.faa'; tbl=pdir/'rvt.tblout'; dom=pdir/'rvt.domtblout'
    run([prodigal,'-i',str(fa),'-o',str(gff),'-a',str(faa),'-f','gff','-p','meta','-q'])
    run([hmmscan,'--noali','--cpu',str(cpu),'-E',str(EMAX),'--domE',str(DOMMAX),'--tblout',str(tbl),'--domtblout',str(dom),str(rvt),str(faa)])
    coords=parse_prodigal_faa_coords(faa)
    for q,model,fe,ie,fs,ds in parse_domtblout(dom):
        if q not in coords:continue
        s,t,strand=coords[q]; nh,rq,rp=local_rec(s,t,rec)
        rows.append({'case_id':cid,'method':'PRODIGAL_RVT','candidate_id':q,'start0':s,'end0':t,'strand':strand,'hmm_model':model,'full_evalue':fe,'domain_i_evalue':ie,'hmm_score':fs,'domain_score':ds,'local_recurrence_hit_n':nh,'best_local_recurrence_q':rq,'recurrence_source':rp})
    sfdir=cdir/'sixframe'; sfdir.mkdir(exist_ok=True); sffaa=sfdir/'sixframe.faa'; sftbl=sfdir/'rvt.tblout'; sfdom=sfdir/'rvt.domtblout'; coord={}
    recs=list(fasta_records(fa)); assert len(recs)==1; seq=recs[0][1]
    with sffaa.open('w') as f:
        for i,(strand,frame,a,b,s,t,pep) in enumerate(sixframe(seq,MINAA),1):
            q=f'SF{i:08d}'; coord[q]=(s,t,strand); f.write(f'>{q} start0={s} end0={t} strand={strand} frame={frame}\n{pep}\n')
    run([hmmscan,'--noali','--cpu',str(cpu),'-E',str(EMAX),'--domE',str(DOMMAX),'--tblout',str(sftbl),'--domtblout',str(sfdom),str(rvt),str(sffaa)])
    for q,model,fe,ie,fs,ds in parse_domtblout(sfdom):
        if q not in coord:continue
        s,t,strand=coord[q]; nh,rq,rp=local_rec(s,t,rec)
        rows.append({'case_id':cid,'method':'SIXFRAME_RVT','candidate_id':q,'start0':s,'end0':t,'strand':strand,'hmm_model':model,'full_evalue':fe,'domain_i_evalue':ie,'hmm_score':fs,'domain_score':ds,'local_recurrence_hit_n':nh,'best_local_recurrence_q':rq,'recurrence_source':rp})
    fields=['case_id','method','candidate_id','start0','end0','strand','hmm_model','full_evalue','domain_i_evalue','hmm_score','domain_score','local_recurrence_hit_n','best_local_recurrence_q','recurrence_source']
    write_tsv(cdir/'BASELINE_RT_CANDIDATES.tsv',rows,fields)
    ans={'case_id':cid,'length_bp':n,'recurrence_proposal_n':len(rec)}
    for method in ['PRODIGAL_RVT','SIXFRAME_RVT']:
        rr=[r for r in rows if r['method']==method]; iv=[]
        for r in rr:iv.append((max(0,int(r['start0'])-RADIUS),min(n,int(r['end0'])+RADIUS)))
        ans[method.lower()+'_n']=len(rr); ans[method.lower()+'_candidate_per_mb']=len(rr)/(n/1e6); ans[method.lower()+'_exposure_fraction']=union_len(iv)/n if n else 0
    return ans

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--cpu',type=int,default=16); a=ap.parse_args(); here=Path(__file__).resolve().parents[1]
    gate=json.loads((here/'frozen/BASELINE_B0_GATE.json').read_text())
    if not gate['gate_pass']:raise SystemExit('B0 R1 gate is HOLD; do not score.')
    cases=read_tsv(here/'config/CASES.tsv'); seq={r['case_id']:r for r in read_tsv(here/'frozen/BASELINE_SEQUENCE_MANIFEST.tsv')}; tool=json.loads((here/'frozen/BASELINE_TOOLCHAIN.json').read_text())
    rvt=Path(tool['rvt_path']); prodigal=tool['tools']['prodigal']['path']; hmmscan=tool['tools']['hmmscan']['path']
    out=here/'outputs/score_freeze'; out.mkdir(parents=True,exist_ok=True); summary=[]
    # CRITICAL: this runner never reads BASELINE_TRUTH_SEAL.tsv or any truth artifact.
    for c in cases:
        z=score_case(c,seq[c['case_id']],rvt,prodigal,hmmscan,out,a.cpu,here); summary.append(z); print('[SCORE]',z,flush=True)
    write_tsv(out/'BASELINE_SCORE_SUMMARY.tsv',summary,list(summary[0]))
    hashes=[]
    for p in sorted(out.rglob('*')):
        if p.is_file():hashes.append({'path':str(p.relative_to(here)),'sha256':sha256_file(p)})
    write_tsv(out/'BASELINE_SCORE_FREEZE_SHA256.tsv',hashes,['path','sha256'])
    gate={'stage':'Stage05F4A-B1-R1-score-freeze','version':'v1.0.1','gate':'PASS_STAGE05F4A_B1_R1_PROTEIN_FIRST_BASELINE_SCORES_FROZEN_NO_TRUTH','gate_pass':True,'truth_opened_for_scoring':False,'retuning_authorized':False,'rt_neighborhood_radius_bp':RADIUS,'rvt_full_evalue':EMAX,'rvt_domain_i_evalue':DOMMAX,'sixframe_min_aa':MINAA,'prodigal_path':prodigal,'hmmscan_path':hmmscan}
    safe_json(out/'BASELINE_B1_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__':main()
