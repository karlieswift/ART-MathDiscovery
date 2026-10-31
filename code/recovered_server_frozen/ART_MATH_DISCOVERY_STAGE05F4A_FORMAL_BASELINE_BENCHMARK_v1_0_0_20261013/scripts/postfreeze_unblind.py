#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math, re
from pathlib import Path
from common import *

def ov(a,b,c,d): return a<d and b>c

def rank_candidates(rows, method):
    x=[r for r in rows if r['method']==method]
    # RT-only rank: HMM evidence; +REC rank: recurrence first then HMM
    rt=sorted(x,key=lambda r:(float(r['hmm_evalue']),-float(r['hmm_score'])))
    rr=sorted(x,key=lambda r:(0 if int(r['local_recurrence_hit_n'])>0 else 1,float(r['best_local_recurrence_q']),float(r['hmm_evalue']),-float(r['hmm_score'])))
    return rt,rr

def first_truth_rank(rows,t0,t1):
    for i,r in enumerate(rows,1):
        if ov(int(r['start0']),int(r['end0']),t0,t1): return i
    return None

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--truth-tsv',default=None); a=ap.parse_args()
    here=Path(__file__).resolve().parents[1]
    b1=json.loads((here/'outputs/score_freeze/BASELINE_B1_GATE.json').read_text())
    if not b1['gate_pass']: raise SystemExit('B1 not frozen')
    truth=Path(a.truth_tsv) if a.truth_tsv else here/'frozen/BASELINE_TRUTH_SEAL.tsv'
    cases=read_tsv(truth); out=[]; unresolved=[]
    for c in cases:
        if c['primary_eval']!='1': continue
        if c['truth_start0'] in {'','.','NA'} or c['truth_end0'] in {'','.','NA'}:
            unresolved.append(c['case_id']); continue
        t0=int(c['truth_start0']); t1=int(c['truth_end0'])
        p=here/'outputs/score_freeze'/c['case_id']/'BASELINE_RT_CANDIDATES.tsv'
        rows=read_tsv(p)
        for method in ['PRODIGAL_RVT','SIXFRAME_RVT']:
            rt,rr=rank_candidates(rows,method)
            for mode,arr in [('RVT_ONLY',rt),('RVT_PLUS_FROZEN_REC',rr)]:
                rank=first_truth_rank(arr,t0,t1); n=len(arr)
                out.append({'case_id':c['case_id'],'group':c['group'],'method':method+'_'+mode,'truth_start0':t0,'truth_end0':t1,'candidate_n':n,'truth_rank':rank or '', 'recovered':int(rank is not None),'hit1':int(rank==1),'hit3':int(rank is not None and rank<=3),'hit5':int(rank is not None and rank<=5),'hit10':int(rank is not None and rank<=10),'reciprocal_rank':(1.0/rank if rank else 0.0),'normalized_rank_percentile':(100.0*rank/n if rank and n else '')})
    fields=list(out[0]) if out else ['case_id']
    write_tsv(here/'outputs/BASELINE_POSTFREEZE_EVALUATION.tsv',out,fields)
    gate_pass=(len(unresolved)==0)
    gate={'stage':'Stage05F4A-B2-postfreeze-evaluation','gate':'PASS_STAGE05F4A_B2_PRIMARY_TRUTH_EVALUATED' if gate_pass else 'HOLD_STAGE05F4A_B2_UNRESOLVED_PROJECT_TRUTH','gate_pass':gate_pass,'unresolved_truth_cases':unresolved,'truth_used_before_score_freeze':False,'retuning_authorized':False}
    safe_json(here/'outputs/BASELINE_B2_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
