#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd

EXPECTED_SOURCE_SHA256 = "272b59fd09d853bc7336e364a64e909c3cca7b082a684c69355d873707be9337"
EXPECTED_UNIVERSE_SHA256 = "24f69c4f463678413d772cc567332bd0b94ff18019d53ea7e454aa5e30fb5ea4"
KEEP = [
    "analysis_id","sequence_id","arm","length_bp","region_rank","region_id","start0","end0","span_bp",
    "selected_window_rows","best_q_acnv_family","best_q_recurrence","sum_family_surprisal","best_single_family_q",
    "acnv_best_member","acnv_best_window_start0","acnv_best_window_end0","acnv_best_window_len",
    "recurrence_best_k","recurrence_best_lag_bp","recurrence_best_window_start0","recurrence_best_window_end0",
    "recurrence_best_window_len"
]

def sha256(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()

def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--source', required=True, type=Path)
    ap.add_argument('--universe', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--allow-noncanonical-source', action='store_true', help='Selftest only; never use for formal closeout.')
    a=ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    src_sha, uni_sha = sha256(a.source), sha256(a.universe)
    if not a.allow_noncanonical_source:
        if src_sha != EXPECTED_SOURCE_SHA256:
            raise SystemExit(f"[FAIL] frozen region board SHA256 mismatch: {src_sha}")
        if uni_sha != EXPECTED_UNIVERSE_SHA256:
            raise SystemExit(f"[FAIL] universe manifest SHA256 mismatch: {uni_sha}")
    df=pd.read_csv(a.source,sep='\t')
    uni=pd.read_csv(a.universe,sep='\t')
    missing=[c for c in KEEP+['support_families'] if c not in df.columns]
    if missing: raise SystemExit(f"[FAIL] missing columns: {missing}")
    s=df['support_families'].fillna('').astype(str)
    dual=df.loc[s.str.contains('ACNV_FAMILY',regex=False) & s.str.contains('GENERIC_RECURRENCE',regex=False), KEEP].copy()
    dual=dual.sort_values(['sequence_id','region_rank','start0','end0']).reset_index(drop=True)
    dual['two_arm_max_q']=dual[['best_q_acnv_family','best_q_recurrence']].max(axis=1)
    dual['two_arm_min_surprisal']=-np.log10(dual['two_arm_max_q'].clip(lower=np.nextafter(0,1)))
    dual.to_csv(a.outdir/'STAGE04A_R1_DUAL_ARM_REGIONS.tsv',sep='\t',index=False,float_format='%.12g')
    per=uni[['analysis_id','sequence_id','arm','length_bp']].copy()
    cnt=dual.groupby('sequence_id').size().rename('dual_arm_region_n')
    br=dual.groupby('sequence_id')['region_rank'].min().rename('best_dual_arm_region_rank')
    per=per.merge(cnt,left_on='sequence_id',right_index=True,how='left').merge(br,left_on='sequence_id',right_index=True,how='left')
    per['dual_arm_region_n']=per['dual_arm_region_n'].fillna(0).astype(int)
    per['best_dual_arm_region_rank']=per['best_dual_arm_region_rank'].astype('Int64')
    per['dual_arm_burden_per_Mb']=per['dual_arm_region_n']/(per['length_bp']/1e6)
    per.to_csv(a.outdir/'STAGE04A_R1_PER_SEQUENCE_DUAL_ARM_BURDEN.tsv',sep='\t',index=False,float_format='%.12g')
    best=(dual.sort_values(['sequence_id','region_rank','two_arm_max_q','start0','end0']).groupby('sequence_id',as_index=False).first())
    best.to_csv(a.outdir/'STAGE04A_R1_BEST_DUAL_ARM_PER_SEQUENCE.tsv',sep='\t',index=False,float_format='%.12g')
    total_mb=float(uni['length_bp'].sum()/1e6)
    n_all=len(df); n_dual=len(dual)
    formal = (src_sha == EXPECTED_SOURCE_SHA256 and uni_sha == EXPECTED_UNIVERSE_SHA256)
    summary={
      'stage':'Stage04A-R1',
      'gate':'PASS_STAGE04A_R1_FROZEN_TWO_ARM_UNIVERSE_CLOSEOUT' if formal else 'PASS_SELFTEST_NONCANONICAL_SOURCE',
      'gate_pass':True,
      'formal_canonical_source':formal,
      'source_frozen_region_board':a.source.name,
      'source_frozen_region_board_sha256':src_sha,
      'source_universe_manifest_sha256':uni_sha,
      'n_sequences':int(len(uni)),
      'n_all_frozen_union_regions':int(n_all),
      'n_dual_arm_regions':int(n_dual),
      'dual_arm_fraction_of_all_regions':float(n_dual/n_all),
      'candidate_space_compression_fold':float(n_all/n_dual),
      'n_sequences_with_ge1_dual_arm':int((per.dual_arm_region_n>0).sum()),
      'fraction_sequences_with_ge1_dual_arm':float((per.dual_arm_region_n>0).mean()),
      'n_sequences_without_dual_arm':int((per.dual_arm_region_n==0).sum()),
      'total_Mb':total_mb,
      'pooled_dual_arm_burden_per_Mb':float(n_dual/total_mb),
      'median_dual_arm_regions_per_sequence':float(per.dual_arm_region_n.median()),
      'iqr_dual_arm_regions_per_sequence':[float(per.dual_arm_region_n.quantile(.25)),float(per.dual_arm_region_n.quantile(.75))],
      'p95_dual_arm_regions_per_sequence':float(per.dual_arm_region_n.quantile(.95)),
      'median_dual_arm_burden_per_Mb':float(per.dual_arm_burden_per_Mb.median()),
      'iqr_dual_arm_burden_per_Mb':[float(per.dual_arm_burden_per_Mb.quantile(.25)),float(per.dual_arm_burden_per_Mb.quantile(.75))],
      'p95_dual_arm_burden_per_Mb':float(per.dual_arm_burden_per_Mb.quantile(.95)),
      'dual_arm_region_rank_median':float(dual.region_rank.median()),
      'dual_arm_region_rank_iqr':[float(dual.region_rank.quantile(.25)),float(dual.region_rank.quantile(.75))],
      'n_sequences_best_dual_rank_eq1':int((per.best_dual_arm_region_rank==1).sum()),
      'n_sequences_best_dual_rank_le2':int((per.best_dual_arm_region_rank<=2).sum()),
      'n_sequences_best_dual_rank_le3':int((per.best_dual_arm_region_rank<=3).sum()),
      'n_sequences_best_dual_rank_le5':int((per.best_dual_arm_region_rank<=5).sum()),
      'n_sequences_best_dual_rank_le10':int((per.best_dual_arm_region_rank<=10).sum()),
      'score_recomputed':False,'region_ranking_recomputed':False,'threshold_changed':False,
      'annotation_used':False,'truth_used':False,'retuning_authorized':False,
      'candidate_interpretation':'Mathematical two-arm candidate regions for downstream biological annotation; not ART identity calls.'
    }
    (a.outdir/'STAGE04A_R1_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    gate={k:summary[k] for k in ['stage','gate','gate_pass','formal_canonical_source','source_frozen_region_board_sha256','source_universe_manifest_sha256','n_sequences','n_all_frozen_union_regions','n_dual_arm_regions','score_recomputed','region_ranking_recomputed','threshold_changed','annotation_used','truth_used','retuning_authorized','candidate_interpretation']}
    (a.outdir/'STAGE04A_R1_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return 0
if __name__=='__main__': raise SystemExit(main())
