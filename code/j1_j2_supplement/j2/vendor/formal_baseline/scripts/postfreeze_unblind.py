#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
from common import *

def ov(a,b,c,d):return a<d and b>c

def truth_rank(rows,t0,t1,start='start0',end='end0'):
    for i,r in enumerate(rows,1):
        if ov(int(r[start]),int(r[end]),t0,t1):return i
    return None

def rt_orders(rows,method):
    x=[r for r in rows if r['method']==method]
    rt=sorted(x,key=lambda r:(float(r['full_evalue']),float(r['domain_i_evalue']),-float(r['hmm_score'])))
    rr=sorted(x,key=lambda r:(0 if int(r['local_recurrence_hit_n'])>0 else 1,float(r['best_local_recurrence_q']),float(r['full_evalue']),float(r['domain_i_evalue']),-float(r['hmm_score'])))
    return rt,rr

def add(out,c,method,arr,t0,t1):
    r=truth_rank(arr,t0,t1); n=len(arr)
    out.append({'case_id':c['case_id'],'group':c['group'],'role':c['role'],'method':method,'candidate_n':n,'truth_rank':r or '','recovered':int(r is not None),'hit1':int(r==1),'hit3':int(r is not None and r<=3),'hit5':int(r is not None and r<=5),'hit10':int(r is not None and r<=10),'reciprocal_rank':1.0/r if r else 0.0,'normalized_rank_percentile':100.0*r/n if r and n else ''})

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--truth-tsv',default=None); a=ap.parse_args(); here=Path(__file__).resolve().parents[1]
    b1=json.loads((here/'outputs/score_freeze/BASELINE_B1_GATE.json').read_text())
    if not b1['gate_pass']:raise SystemExit('B1 not frozen')
    truth=Path(a.truth_tsv) if a.truth_tsv else here/'frozen/BASELINE_TRUTH_SEAL.tsv'; cases=read_tsv(truth); out=[]; unresolved=[]
    for c in cases:
        if c['primary_eval']!='1':continue
        if c['truth_start0'] in {'','.','NA'} or c['truth_end0'] in {'','.','NA'}:
            unresolved.append(c['case_id']);continue
        t0=int(c['truth_start0']); t1=int(c['truth_end0']); cid=c['case_id']
        # Frozen raw-DNA methods: canonicalized before truth opening.
        props=read_tsv(here/'frozen/canonical'/cid/'CANONICAL_PROPOSALS.tsv'); unions=read_tsv(here/'frozen/canonical'/cid/'CANONICAL_UNION_REGIONS.tsv')
        unions=sorted(unions,key=lambda r:(int(r['source_rank']),int(r['source_row'])))
        ac=sorted([r for r in props if 'ACNV' in r['family']],key=lambda r:(float(r['qvalue']),int(r['source_row'])))
        rr=sorted([r for r in props if 'RECURRENCE' in r['family']],key=lambda r:(float(r['qvalue']),int(r['source_row'])))
        add(out,c,'ARM_M_FULL_FROZEN',unions,t0,t1); add(out,c,'ACNV_ONLY_FROZEN',ac,t0,t1); add(out,c,'RECURRENCE_ONLY_FROZEN',rr,t0,t1)
        rows=read_tsv(here/'outputs/score_freeze'/cid/'BASELINE_RT_LOCI.tsv')
        for method in ['PRODIGAL_RVT','SIXFRAME_RVT']:
            rt,rrec=rt_orders(rows,method); add(out,c,method+'_ONLY',rt,t0,t1); add(out,c,method+'_PLUS_FROZEN_REC',rrec,t0,t1)
    fields=list(out[0]) if out else ['case_id']; write_tsv(here/'outputs/BASELINE_POSTFREEZE_EVALUATION.tsv',out,fields)
    gate_pass=not unresolved
    gate={'stage':'Stage05F4A-B2-R2-postfreeze-evaluation','version':'v1.0.2','gate':'PASS_STAGE05F4A_B2_R2_PRIMARY_TRUTH_EVALUATED' if gate_pass else 'HOLD_STAGE05F4A_B2_R2_UNRESOLVED_PROJECT_TRUTH','gate_pass':gate_pass,'unresolved_truth_cases':unresolved,'truth_used_before_score_freeze':False,'retuning_authorized':False}
    safe_json(here/'outputs/BASELINE_B2_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__':main()
