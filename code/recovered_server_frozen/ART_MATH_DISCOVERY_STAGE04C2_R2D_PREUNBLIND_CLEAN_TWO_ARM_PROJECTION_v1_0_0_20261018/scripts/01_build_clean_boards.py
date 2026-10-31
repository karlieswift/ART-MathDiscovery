#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, sys
from pathlib import Path
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument('--root',required=True)
ap.add_argument('--r2c',required=True)
ap.add_argument('--historical-clean',required=True)
a=ap.parse_args()
root=Path(a.root).resolve(); r2c=Path(a.r2c).resolve(); hist=Path(a.historical_clean).resolve()
out=root/'outputs'; croot=out/'contigs'; croot.mkdir(parents=True,exist_ok=True)

# Import exact historical common.py from the frozen historical clean replay root.
spec=importlib.util.spec_from_file_location('historical_clean_common',hist/'scripts/common.py')
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
build_two_arm_board=mod.build_two_arm_board
sha256_file=mod.sha256_file

pol=json.loads((hist/'config/TWO_ARM_POLICY.json').read_text())
thr=float(pol['proposal_q_threshold'])
if thr!=0.005: raise SystemExit('HOLD: historical threshold changed')

lock=list(csv.DictReader((r2c/'outputs/CONTIG_INPUT_LOCK.tsv').open(),delimiter='\t'))
all_regs=[]; all_dual=[]; gates=[]

for x in lock:
    aid=x['assembly_id']; sid=x['scoring_id']
    src=r2c/'outputs/contigs'/sid/'STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz'
    d=pd.read_csv(src,sep='\t')
    if 'sequence_id' not in d.columns: raise SystemExit(f'HOLD missing sequence_id in {sid}')
    seq_ids=set(d['sequence_id'].astype(str).unique())
    if seq_ids!={sid}: raise SystemExit(f'HOLD unexpected sequence_id values in {sid}: {seq_ids}')
    sel,regs=build_two_arm_board(d,thr=thr,prefix='AR')
    cd=croot/sid; cd.mkdir(parents=True,exist_ok=True)
    sp=cd/'CLEAN_TWO_ARM_SELECTED_PROPOSALS.tsv'
    rp=cd/'CLEAN_TWO_ARM_UNION_REGIONS.tsv'
    dp=cd/'CLEAN_TWO_ARM_DUAL_ARM_REGIONS.tsv'
    sel.to_csv(sp,sep='\t',index=False)
    regs.to_csv(rp,sep='\t',index=False)
    dual=regs.loc[regs['support_family_n'].astype(int)==2].copy() if len(regs) else regs.copy()
    dual.to_csv(dp,sep='\t',index=False)
    reg2=regs.copy(); reg2.insert(0,'scoring_id',sid); reg2.insert(0,'assembly_id',aid); all_regs.append(reg2)
    dual2=dual.copy(); dual2.insert(0,'scoring_id',sid); dual2.insert(0,'assembly_id',aid); all_dual.append(dual2)
    gate={
      'gate':'PASS_CLEAN_TWO_ARM_BOARD_FROZEN_PREUNBLIND', 'gate_pass':True,
      'assembly_id':aid,'scoring_id':sid,'proposal_q_threshold':thr,
      'families':['ACNV_FAMILY','GENERIC_RECURRENCE'],
      'n_genomic_windows':int(len(d)),'n_selected_proposals':int(len(sel)),
      'n_union_regions':int(len(regs)),'n_dual_arm_regions':int(len(dual)),
      'source_frozen_scores':str(src),'source_frozen_scores_sha256':sha256_file(src),
      'selected_sha256':sha256_file(sp),'union_regions_sha256':sha256_file(rp),'dual_arm_sha256':sha256_file(dp),
      'truth_used':False,'retuning_authorized':False,'recompute_sequence_features':False,
      'score_recomputed':False,'q_recomputed':False,'science_change':False,
    }
    (cd/'CLEAN_TWO_ARM_BOARD_GATE.json').write_text(json.dumps(gate,indent=2,sort_keys=True)+'\n')
    gates.append(gate)
    print(json.dumps(gate,sort_keys=True))

A=pd.concat(all_regs,ignore_index=True) if all_regs else pd.DataFrame()
D=pd.concat(all_dual,ignore_index=True) if all_dual else pd.DataFrame()
apath=out/'CLEAN_TWO_ARM_ALL_REGIONS.tsv'; dpath=out/'CLEAN_TWO_ARM_DUAL_ARM_VIEW.tsv'
A.to_csv(apath,sep='\t',index=False); D.to_csv(dpath,sep='\t',index=False)
manifest={
 'gate':'PASS_ALL_15_PREUNBLIND_CLEAN_TWO_ARM_BOARDS_FROZEN','gate_pass':True,
 'n_contigs':len(gates),'n_all_clean_regions':int(len(A)),'n_clean_dual_arm_regions':int(len(D)),
 'cases':gates,'truth_used':False,'retuning_authorized':False,'science_change':False,
 'all_regions_sha256':sha256_file(apath),'dual_arm_view_sha256':sha256_file(dpath),
}
(out/'CLEAN_TWO_ARM_BOARD_MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
print(json.dumps(manifest,indent=2,sort_keys=True))
