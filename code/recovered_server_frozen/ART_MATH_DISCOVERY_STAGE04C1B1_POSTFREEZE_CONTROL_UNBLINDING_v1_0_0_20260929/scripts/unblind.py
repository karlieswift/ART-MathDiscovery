#!/usr/bin/env python3
import argparse,csv,gzip,json,hashlib,math
from pathlib import Path

def readj(p): return json.loads(Path(p).read_text())
def readtsv(p):
    with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def ov(a,b,c,d): return int(a)<int(d) and int(c)<int(b)
def f(x):
    try:return float(x)
    except:return math.nan
def tied_rank(vals,target):
    vals=[x for x in vals if math.isfinite(x)]
    lt=sum(x<target for x in vals); le=sum(x<=target for x in vals)
    return lt+1,le,len(vals)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def score_control(c1,cid,truth):
    u=readtsv(c1/'outputs'/cid/'C1_UNION_REGIONS.tsv')
    scorep=c1/'outputs'/cid/'C1_WINDOW_FAMILY_SCORES.tsv.gz'
    rows=[]
    with gzip.open(scorep,'rt',newline='') as g: rows=list(csv.DictReader(g,delimiter='\t'))
    s0,e0=int(truth['start0']),int(truth['end0'])
    ur=[r for r in u if ov(r['start0'],r['end0'],s0,e0)]
    ur=sorted(ur,key=lambda r:int(r['region_rank']))
    ac=[r for r in rows if ov(r['start0'],r['end0'],s0,e0) and math.isfinite(f(r['q_acnv_family']))]
    re=[r for r in rows if ov(r['start0'],r['end0'],s0,e0) and math.isfinite(f(r['q_recurrence']))]
    all_ac=[f(r['q_acnv_family']) for r in rows]; all_re=[f(r['q_recurrence']) for r in rows]
    out={'control_id':cid,'truth_component':truth['component'],'truth_start0':s0,'truth_end0':e0,'truth_span_bp':e0-s0,'n_union_regions':len(u),'truth_overlapping_union_region_n':len(ur)}
    if ur:
        b=ur[0]; out.update(region_rank=int(b['region_rank']),region_id=b['region_id'],region_start0=int(b['start0']),region_end0=int(b['end0']),region_support_family_n=int(b['support_family_n']),region_support_families=b['support_families'],region_top_fraction=int(b['region_rank'])/len(u))
    else: out.update(region_rank='',region_id='',region_start0='',region_end0='',region_support_family_n=0,region_support_families='',region_top_fraction='')
    if ac:
        b=min(ac,key=lambda r:f(r['q_acnv_family'])); q=f(b['q_acnv_family']); rmin,rmax,n=tied_rank(all_ac,q)
        out.update(acnv_best_q=q,acnv_rank_min=rmin,acnv_rank_max=rmax,acnv_rank_den=n,acnv_window_start0=int(b['start0']),acnv_window_end0=int(b['end0']),acnv_best_member=b['acnv_family_best_member'])
    else: out.update(acnv_best_q='',acnv_rank_min='',acnv_rank_max='',acnv_rank_den='',acnv_window_start0='',acnv_window_end0='',acnv_best_member='')
    if re:
        b=min(re,key=lambda r:f(r['q_recurrence'])); q=f(b['q_recurrence']); rmin,rmax,n=tied_rank(all_re,q)
        out.update(recurrence_best_q=q,recurrence_rank_min=rmin,recurrence_rank_max=rmax,recurrence_rank_den=n,recurrence_window_start0=int(b['start0']),recurrence_window_end0=int(b['end0']),recurrence_best_k=int(b['recurrence_best_k']),recurrence_best_lag_bp=int(b['recurrence_best_lag_bp']),recurrence_best_match_count=int(b['recurrence_best_match_count']))
    else: out.update(recurrence_best_q='',recurrence_rank_min='',recurrence_rank_max='',recurrence_rank_den='',recurrence_window_start0='',recurrence_window_end0='',recurrence_best_k='',recurrence_best_lag_bp='',recurrence_best_match_count='')
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--c1',required=True);ap.add_argument('--c1b0',required=True);a=ap.parse_args()
    R=Path(a.root);c1=Path(a.c1);b0=Path(a.c1b0)
    g1=readj(c1/'outputs/C1_FINAL_GATE.json'); g0=readj(b0/'outputs/C1B0_FINAL_GATE.json')
    if g1.get('gate')!='PASS_STAGE04C1_RT_CONTROL_TWO_ARM_SCORE_FREEZE_COMPLETE':raise RuntimeError('C1 gate mismatch')
    if g0.get('gate')!='PASS_STAGE04C1B0_RT_ANCHOR_TRUTH_RESOLUTION':raise RuntimeError('C1B0 gate mismatch')
    truth=readtsv(b0/'outputs/C1B0_CONTROL_TRUTH_RESOLVED.tsv')
    prim=[r for r in truth if r['component'].startswith('RT_')]
    if len(prim)!=3:raise RuntimeError(f'expected 3 primary RT anchors, got {len(prim)}')
    outs=[score_control(c1,r['control_id'],r) for r in prim]
    cols=list(outs[0].keys());op=R/'outputs/C1B1_PRIMARY_RT_ANCHOR_RESULTS.tsv';op.parent.mkdir(exist_ok=True)
    with op.open('w',newline='') as f0:
        w=csv.DictWriter(f0,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(outs)
    # secondary system envelopes
    sec=[]
    for r in truth:
        if r['component']=='CODING_SYSTEM_ENVELOPE':sec.append(score_control(c1,r['control_id'],r))
    if sec:
        cp=R/'outputs/C1B1_SECONDARY_CODING_ENVELOPE_RESULTS.tsv'
        with cp.open('w',newline='') as f0:
            w=csv.DictWriter(f0,fieldnames=list(sec[0].keys()),delimiter='\t');w.writeheader();w.writerows(sec)
    summary=[]
    for x in outs:
        summary.append({'control_id':x['control_id'],'region_rank':x['region_rank'],'n_union_regions':x['n_union_regions'],'dual_arm_at_truth_region':bool(x['region_support_family_n']==2),'acnv_rank':f"{x['acnv_rank_min']}-{x['acnv_rank_max']}/{x['acnv_rank_den']}" if x['acnv_rank_min']!='' else None,'recurrence_rank':f"{x['recurrence_rank_min']}-{x['recurrence_rank_max']}/{x['recurrence_rank_den']}" if x['recurrence_rank_min']!='' else None,'recurrence_best_lag_bp':x['recurrence_best_lag_bp']})
    gate={'gate':'PASS_STAGE04C1B1_POSTFREEZE_CONTROL_UNBLINDING','gate_pass':True,'scientific_status':'POST_SCORE_FREEZE_MECHANISTIC_CONTROL_UNBLINDING','controls_unblinded_n':len(outs),'score_recomputed':False,'region_ranking_recomputed':False,'retuning_authorized':False,'primary_results_sha256':sha(op),'summary':summary}
    (R/'outputs/C1B1_FINAL_GATE.json').write_text(json.dumps(gate,indent=2))
    print(json.dumps(gate,indent=2))
if __name__=='__main__':main()
