#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, json, os, shutil, subprocess, sys, tempfile
from common import read_tsv, write_tsv, dump_json, parse_fasta, revcomp, normalize_system_label, sha256

SYSTEMS=['Alistipes-putredinis','Bacteroides-coprophilus','Bacteroides-sp-9-1-42FAA','Bacteroides-sp-2-1-16','Blautia-hansenii']

def load_cassettes(r3root):
    p=Path(r3root)/'outputs/UG27_R3_CASSETTE_TRUTH_LOCK.tsv'
    rows=read_tsv(p)
    out={}
    for r in rows:
        system=normalize_system_label(r['system'])
        out[system]=r
    return out

def target_fna_for_system(ugroot, system):
    mapping={
      'Alistipes-putredinis':'Alistipes-putredinis-DSM-17216-MAF-2.fna',
      'Bacteroides-coprophilus':'Bacteroides-coprophilus-DSM-18228-MAF-2.fna',
      'Bacteroides-sp-9-1-42FAA':'Bacteroides-sp-9-1-42FAA-MAF-2.fna',
      'Bacteroides-sp-2-1-16':'Bacteroides-sp-2-1-16-MAF-2.fna',
      'Blautia-hansenii':'Blautia-hansenii-DSM-20583-MAF-2.fna'
    }
    return Path(ugroot)/'work/hcom2_7734779/extracted/hCom2/20221117/genomes'/mapping[system]

def get_contig_seq(fna, contig):
    for n,s in parse_fasta(fna):
        if n==contig: return s
    raise KeyError(f'{contig} not found in {fna}')

def resolve_from_tsv(path,cassettes):
    rows=read_tsv(path); req={'system','source_contig_id','start0','end0','strand','unit_id'}
    if not rows or not req.issubset(rows[0]): raise ValueError('truth TSV needs columns '+','.join(sorted(req)))
    out=[]
    for r in rows:
        s=normalize_system_label(r['system'])
        if s not in cassettes: continue
        c=cassettes[s]; a0=int(r['start0']); a1=int(r['end0'])
        cs=int(c['cassette_start0']); ce=int(c['cassette_end0'])
        if r['source_contig_id'] != c['source_contig_id']: raise ValueError(f'contig mismatch for {s}')
        if not (cs <= a0 < a1 <= ce): raise ValueError(f'unit outside frozen cassette for {s}: {a0}-{a1} vs {cs}-{ce}')
        o=dict(r); o['system']=s; o['start0']=a0; o['end0']=a1; o['source_method']=r.get('source_method','AUTHOR_MACHINE_READABLE_EXACT_UNIT_COORDINATES'); o['evalue']=r.get('evalue','')
        out.append(o)
    return out,'AUTHOR_MACHINE_READABLE_EXACT_UNIT_COORDINATES'

def parse_tblout(path, model_name, system, contig):
    hits=[]
    with open(path,encoding='utf-8',errors='ignore') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            x=line.split()
            # Infernal --tblout: target name acc query name acc mdl mdl_from mdl_to seq_from seq_to strand trunc pass gc bias score E inc description
            if len(x) < 18: continue
            target=x[0]; query=x[2]
            try:
                sf=int(x[7]); st=int(x[8]); strand=x[9]; score=float(x[14]); ev=float(x[15])
            except Exception: continue
            start0=min(sf,st)-1; end0=max(sf,st)
            hits.append({'system':system,'source_contig_id':contig,'start0':start0,'end0':end0,'strand':strand,'unit_id':'','evalue':ev,'cm_score':score,'model':model_name,'query_model':query,'source_method':'AUTHOR_FINAL_UG27_COVARIANCE_MODELS_PUBLISHED_CMSEARCH_RULE'})
    return hits

def dedup_model_overlap(hits):
    hits=sorted(hits,key=lambda h:(h['system'],h['source_contig_id'],h['start0'],h['end0'],h['evalue']))
    keep=[]
    for h in hits:
        conflict=[]
        for i,k in enumerate(keep):
            if h['system']==k['system'] and h['source_contig_id']==k['source_contig_id'] and max(h['start0'],k['start0']) < min(h['end0'],k['end0']):
                conflict.append(i)
        if not conflict:
            keep.append(h); continue
        candidates=[h]+[keep[i] for i in conflict]
        best=min(candidates,key=lambda z:(float(z['evalue']),-float(z.get('cm_score',0)),z['start0']))
        for i in sorted(conflict,reverse=True): keep.pop(i)
        keep.append(best)
    keep=sorted(keep,key=lambda h:(h['system'],h['start0'],h['end0']))
    counters={}
    for h in keep:
        counters[h['system']]=counters.get(h['system'],0)+1
        h['unit_id']=f"{h['system']}_U{counters[h['system']]:02d}"
    return keep

def resolve_from_cms(cm_paths, ugroot, cassettes, work):
    cmsearch=shutil.which('cmsearch')
    if not cmsearch: raise RuntimeError('cmsearch not found; activate/install an Infernal environment')
    allhits=[]
    for system,c in cassettes.items():
        fna=target_fna_for_system(ugroot,system); contig=c['source_contig_id']; seq=get_contig_seq(fna,contig)
        fa=work/f'{system}.target.fa'; fa.write_text(f'>{contig}\n{seq}\n')
        for cm in cm_paths:
            tbl=work/f'{system}.{Path(cm).stem}.tblout'
            cmd=[cmsearch,'--tblout',str(tbl),'--noali',str(cm),str(fa)]
            p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            (work/f'{system}.{Path(cm).stem}.cmsearch.log').write_text(p.stdout)
            if p.returncode!=0: raise RuntimeError(f'cmsearch failed for {system} {cm}')
            allhits.extend(parse_tblout(tbl,Path(cm).name,system,contig))
    # Published threshold and cassette association.
    filt=[]
    for h in allhits:
        if float(h['evalue']) > 0.01: continue
        c=cassettes[h['system']]; cs=int(c['cassette_start0']); ce=int(c['cassette_end0'])
        if cs <= h['start0'] < h['end0'] <= ce: filt.append(h)
    return dedup_model_overlap(filt),'AUTHOR_FINAL_UG27_COVARIANCE_MODELS_PUBLISHED_CMSEARCH_RULE'

def resolve_from_sequences(path,ugroot,cassettes):
    seqs=list(parse_fasta(path)); out=[]
    for header,useq in seqs:
        # header must begin system|unit_id or system__unit_id
        system=None; unit=header
        for s in SYSTEMS:
            if s in header: system=s; break
        if not system or system not in cassettes: continue
        c=cassettes[system]; contig=c['source_contig_id']; g=get_contig_seq(target_fna_for_system(ugroot,system),contig)
        cs=int(c['cassette_start0']); ce=int(c['cassette_end0']); sub=g[cs:ce]
        hits=[]
        for orient,q in [('+',useq),('-',revcomp(useq))]:
            pos=0
            while True:
                j=sub.find(q,pos)
                if j<0: break
                hits.append((orient,cs+j,cs+j+len(q))); pos=j+1
        if len(hits)!=1: raise ValueError(f'exact unit sequence {header} maps {len(hits)} times within cassette; require unique mapping')
        strand,a0,a1=hits[0]
        out.append({'system':system,'source_contig_id':contig,'start0':a0,'end0':a1,'strand':strand,'unit_id':unit,'evalue':'','model':'','cm_score':'','source_method':'AUTHOR_EXACT_UNIT_SEQUENCES_UNIQUELY_MAPPED_TO_FROZEN_NATIVE_CASSETTE'})
    return out,'AUTHOR_EXACT_UNIT_SEQUENCES_UNIQUELY_MAPPED_TO_FROZEN_NATIVE_CASSETTE'

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ug-root',required=True); ap.add_argument('--r3-root',required=True); ap.add_argument('--out',required=True)
    ap.add_argument('--truth-tsv',default=''); ap.add_argument('--cm-main',default=''); ap.add_argument('--cm-alt',default=''); ap.add_argument('--unit-fasta',default='')
    a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True); work=out/'work_cmsearch'; work.mkdir(exist_ok=True)
    c=load_cassettes(a.r3_root)
    mode=''; units=[]; problems=[]
    try:
        if a.truth_tsv:
            units,mode=resolve_from_tsv(a.truth_tsv,c)
        elif a.cm_main and a.cm_alt:
            units,mode=resolve_from_cms([Path(a.cm_main),Path(a.cm_alt)],a.ug_root,c,work)
        elif a.unit_fasta:
            units,mode=resolve_from_sequences(a.unit_fasta,a.ug_root,c)
        else:
            problems.append('No acceptable machine-readable author exact-truth source supplied. Need exact coordinate table, both final UG27 covariance models, or exact author unit sequences.')
    except Exception as e:
        problems.append(str(e))
    by={s:[] for s in SYSTEMS}
    for u in units:
        if u['system'] in by: by[u['system']].append(u)
    for s in SYSTEMS:
        by[s]=sorted(by[s],key=lambda z:(int(z['start0']),int(z['end0'])))
        # Require at least 2 author-defined units to call an array.
        if mode and len(by[s])<2: problems.append(f'{s}: fewer than 2 exact author-defined ncRNA units ({len(by[s])})')
    fields=['system','source_contig_id','start0','end0','strand','unit_id','evalue','cm_score','model','source_method']
    write_tsv(out/'UG27_EXACT_NCRNA_UNITS.tsv',units,fields)
    spans=[]
    if mode and not problems:
        for s in SYSTEMS:
            uu=by[s]; c0=c[s]
            spans.append({'system':s,'source_contig_id':c0['source_contig_id'],'array_start0':min(int(x['start0']) for x in uu),'array_end0':max(int(x['end0']) for x in uu),'array_span_bp':max(int(x['end0']) for x in uu)-min(int(x['start0']) for x in uu),'exact_unit_n':len(uu),'cassette_start0':c0['cassette_start0'],'cassette_end0':c0['cassette_end0'],'truth_source_method':mode})
    write_tsv(out/'UG27_EXACT_NCRNA_ARRAY_SPANS.tsv',spans,['system','source_contig_id','array_start0','array_end0','array_span_bp','exact_unit_n','cassette_start0','cassette_end0','truth_source_method'])
    ok=(mode!='' and not problems and len(spans)==5)
    gate={'stage':'Stage04C2-R4A','gate':'PASS_STAGE04C2_R4A_FIVE_EXACT_AUTHOR_NCRNA_ARRAY_TRUTHS_FROZEN' if ok else 'HOLD_STAGE04C2_R4A_EXACT_ARRAY_TRUTH_NOT_RESOLVED','gate_pass':ok,'truth_source_method':mode or None,'n_prespecified_systems':5,'n_resolved_systems':len(spans),'n_exact_units_total':len(units),'published_cmsearch_unit_evalue_max':0.01,'nearby_threshold_defined':False,'detector_scores_used_for_truth_resolution':False,'detector_regions_used_for_truth_resolution':False,'motif_only_heuristic_used':False,'generic_intergenic_boundary_used_as_exact_truth':False,'retuning_authorized':False,'score_recomputed':False,'science_change':False,'problems':problems}
    if (out/'UG27_EXACT_NCRNA_UNITS.tsv').exists(): gate['exact_units_sha256']=sha256(out/'UG27_EXACT_NCRNA_UNITS.tsv')
    if (out/'UG27_EXACT_NCRNA_ARRAY_SPANS.tsv').exists(): gate['exact_array_spans_sha256']=sha256(out/'UG27_EXACT_NCRNA_ARRAY_SPANS.tsv')
    dump_json(out/'STAGE04C2_R4A_TRUTH_GATE.json',gate); print(json.dumps(gate,indent=2))
    if not ok: sys.exit(3)
if __name__=='__main__': main()
