#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import pandas as pd
from common import sha256_file, strict_overlap, overlap_bp, interval_gap_bp, classify_direct

EXPECTED_ALL='1d8f394a5e35e48301c4c0ebe959b2babb301269eb5977bb9c66c9c86545e3c9'
EXPECTED_DUAL='12940178a8efbff191c066e9bec28ed60ce1b4dad1da7507ba9049995c9dbe15'

def fail(x): raise SystemExit('HOLD_STAGE04C2_R3_UNBLIND: '+x)

def choose_nearest(df,c0,c1):
    if len(df)==0: return None
    z=df.copy()
    z['_gap']=[interval_gap_bp(c0,c1,a,b) for a,b in zip(z.start0,z.end0)]
    z['_ov']=[strict_overlap(c0,c1,a,b) for a,b in zip(z.start0,z.end0)]
    z=z.sort_values(['_gap','region_rank','start0','end0'],ascending=[True,True,True,True])
    return z.iloc[0]

def val(row,key,default=''):
    if row is None or key not in row.index: return default
    x=row[key]
    if pd.isna(x): return default
    return x

def region_payload(prefix,row,c0,c1):
    keys=['region_rank','region_id','start0','end0','span_bp','support_family_n','support_families',
          'best_q_acnv_family','best_q_recurrence','sum_two_family_surprisal','best_single_family_q',
          'acnv_best_member','acnv_best_window_start0','acnv_best_window_end0','acnv_best_window_len',
          'recurrence_best_k','recurrence_best_lag_bp','recurrence_best_window_start0','recurrence_best_window_end0','recurrence_best_window_len']
    out={}
    for k in keys: out[prefix+k]=val(row,k,'')
    if row is None:
        out[prefix+'strict_overlap']=''; out[prefix+'overlap_bp']=''; out[prefix+'gap_bp']=''
    else:
        out[prefix+'strict_overlap']=bool(strict_overlap(c0,c1,int(row.start0),int(row.end0)))
        out[prefix+'overlap_bp']=int(overlap_bp(c0,c1,int(row.start0),int(row.end0)))
        out[prefix+'gap_bp']=int(interval_gap_bp(c0,c1,int(row.start0),int(row.end0)))
    return out

ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--r2d',required=True); a=ap.parse_args()
root=Path(a.root).resolve(); r2d=Path(a.r2d).resolve(); out=root/'outputs'
for p in ['STAGE04C2_R3_SCORE_ARCHIVE_PREFLIGHT_GATE.json','STAGE04C2_R3_TRUTH_LOCK_GATE.json']:
    G=json.loads((out/p).read_text())
    if not G.get('gate_pass'): fail(f'prior gate failed: {p}')
allp=r2d/'outputs/CLEAN_TWO_ARM_ALL_REGIONS.tsv'; dualp=r2d/'outputs/CLEAN_TWO_ARM_DUAL_ARM_VIEW.tsv'
if sha256_file(allp)!=EXPECTED_ALL or sha256_file(dualp)!=EXPECTED_DUAL: fail('R2D clean board hashes changed')
A=pd.read_csv(allp,sep='\t'); D=pd.read_csv(dualp,sep='\t'); T=pd.read_csv(out/'UG27_R3_CASSETTE_TRUTH_LOCK.tsv',sep='\t')
if len(T)!=5: fail('truth lock is not five rows')
summary=[]; details=[]; nearest=[]
for _,t in T.iterrows():
    sid=str(t.scoring_id); c0=int(t.cassette_start0); c1=int(t.cassette_end0)
    a=A.loc[A.scoring_id.astype(str)==sid].copy(); d=D.loc[D.scoring_id.astype(str)==sid].copy()
    if len(a)==0: fail(f'no clean all-region board for truth scoring contig {sid}')
    direct=a.loc[(a.start0.astype(int)<c1)&(a.end0.astype(int)>c0)].copy()
    direct_dual=direct.loc[direct.support_family_n.astype(int)==2].copy()
    direct_acnv=direct.loc[(direct.support_family_n.astype(int)==1)&(direct.support_families.astype(str)=='ACNV_FAMILY')].copy()
    direct_rec=direct.loc[(direct.support_family_n.astype(int)==1)&(direct.support_families.astype(str)=='GENERIC_RECURRENCE')].copy()
    cls=classify_direct(direct)
    if cls=='DUAL_ARM_DIRECT' and len(direct_dual)==0: fail('classification inconsistency')
    primary=bool(len(direct_dual)>0)
    best_direct=direct.sort_values(['region_rank','start0']).iloc[0] if len(direct) else None
    best_direct_dual=direct_dual.sort_values(['region_rank','start0']).iloc[0] if len(direct_dual) else None
    near_all=choose_nearest(a,c0,c1); near_dual=choose_nearest(d,c0,c1)
    for _,r in direct.sort_values(['region_rank','start0']).iterrows():
        details.append({
          'system':t.system,'assembly_id':t.assembly_id,'scoring_id':sid,'cassette_start0':c0,'cassette_end0':c1,
          'cassette_span_bp':int(t.cassette_span_bp),'region_rank':int(r.region_rank),'region_id':r.region_id,
          'region_start0':int(r.start0),'region_end0':int(r.end0),'support_family_n':int(r.support_family_n),
          'support_families':r.support_families,'overlap_bp':int(overlap_bp(c0,c1,int(r.start0),int(r.end0))),
          'best_q_acnv_family':r.best_q_acnv_family,'best_q_recurrence':r.best_q_recurrence,
          'acnv_best_member':r.acnv_best_member,'recurrence_best_k':r.recurrence_best_k,'recurrence_best_lag_bp':r.recurrence_best_lag_bp,
        })
    for kind,r in [('nearest_all_clean',near_all),('nearest_dual_arm',near_dual)]:
        if r is not None:
            nearest.append({
              'system':t.system,'assembly_id':t.assembly_id,'scoring_id':sid,'kind':kind,
              'cassette_start0':c0,'cassette_end0':c1,'region_rank':int(r.region_rank),'region_id':r.region_id,
              'region_start0':int(r.start0),'region_end0':int(r.end0),'support_families':r.support_families,
              'strict_overlap':bool(strict_overlap(c0,c1,int(r.start0),int(r.end0))),
              'overlap_bp':int(overlap_bp(c0,c1,int(r.start0),int(r.end0))),
              'gap_bp':int(interval_gap_bp(c0,c1,int(r.start0),int(r.end0))),
              'best_q_acnv_family':r.best_q_acnv_family,'best_q_recurrence':r.best_q_recurrence,
              'acnv_best_member':r.acnv_best_member,'recurrence_best_k':r.recurrence_best_k,'recurrence_best_lag_bp':r.recurrence_best_lag_bp,
            })
    row={
      'system':t.system,'assembly_id':t.assembly_id,'scoring_id':sid,'source_contig_id':t.source_contig_id,
      'orientation':t.orientation,'pident':float(t.pident),'cassette_start0':c0,'cassette_end0':c1,'cassette_span_bp':int(t.cassette_span_bp),
      'clean_all_region_n_on_contig':int(len(a)),'clean_dual_arm_region_n_on_contig':int(len(d)),
      'direct_clean_region_n':int(len(direct)),'direct_dual_arm_region_n':int(len(direct_dual)),
      'direct_acnv_only_region_n':int(len(direct_acnv)),'direct_recurrence_only_region_n':int(len(direct_rec)),
      'primary_cassette_dual_arm_strict_overlap':primary,'cassette_level_classification':cls,
      'exact_ncrna_array_truth':False,
    }
    row.update(region_payload('best_direct_',best_direct,c0,c1))
    row.update(region_payload('best_direct_dual_',best_direct_dual,c0,c1))
    row.update(region_payload('nearest_all_',near_all,c0,c1))
    row.update(region_payload('nearest_dual_',near_dual,c0,c1))
    summary.append(row)
S=pd.DataFrame(summary).sort_values(['assembly_id','system']).reset_index(drop=True)
O=pd.DataFrame(details)
N=pd.DataFrame(nearest)
sp=out/'UG27_CASSETTE_UNBLINDING_SUMMARY.tsv'; op=out/'UG27_CASSETTE_DIRECT_OVERLAPS.tsv'; np=out/'UG27_CASSETTE_NEAREST_REGIONS.tsv'
S.to_csv(sp,sep='\t',index=False); O.to_csv(op,sep='\t',index=False); N.to_csv(np,sep='\t',index=False)
print(S[['system','assembly_id','scoring_id','cassette_level_classification','primary_cassette_dual_arm_strict_overlap','direct_clean_region_n','direct_dual_arm_region_n','nearest_dual_gap_bp']].to_string(index=False))
