#!/usr/bin/env python3
import argparse,csv,json,tarfile
from pathlib import Path
import pandas as pd
from common import *
ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--marshill-root',required=True);ap.add_argument('--stage03a-tar',required=True);a=ap.parse_args()
root=Path(a.root).resolve();mr=Path(a.marshill_root).resolve();tar=Path(a.stage03a_tar).resolve();boards=root/'outputs/boards'
# enforce all boards existed before truth read
bm=read_json(root/'outputs/CLEAN_TWO_ARM_BOARD_MANIFEST.json')
if not bm.get('gate_pass') or bm.get('truth_used') is not False: raise SystemExit('FAIL boards not cleanly frozen')
truth=[]
mj=read_json(root/'truth/MARSHILL_POSTHOC_TRUTH.json');truth.append({'case':'MW248466.1_MarsHill','accession':'MW248466.1','name':'MarsHill','start0':int(mj['start0']),'end0':int(mj['end0']),'truth_type':mj['truth_interval_type']})
with open(root/'truth/ART_TRUTH_MANIFEST_STAGE03B.tsv') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        truth.append({'case':r['accession'],'accession':r['accession'],'name':r['name'],'start0':int(r['strict_repeat_core_start0']),'end0':int(r['strict_repeat_core_end0']),'truth_type':'Stage03B_strict_repeat_core'})
# frozen score accessors
score_map={'MW248466.1_MarsHill':pd.read_csv(mr/'outputs/STAGE02B_WINDOW_FAMILY_SCORES.tsv.gz',sep='\t')}
with tarfile.open(tar,'r:gz') as tf:
    for acc in ['MW349129','MW349128','MW218148.1']:
        score_map[acc],_=read_tsv_from_tar(tf,f'outputs/{acc}/STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz')
rows=[]
for t in truth:
    regs=pd.read_csv(boards/t['case']/'CLEAN_TWO_ARM_UNION_REGIONS.tsv',sep='\t'); s0,e0=t['start0'],t['end0']
    hits=regs[(regs.start0<e0)&(regs.end0>s0)]
    d=score_map[t['case']]
    ba=best_truth_window(d,s0,e0,'q_acnv_family',[('best_member','acnv_family_best_member','str')])
    br=best_truth_window(d,s0,e0,'q_recurrence',[('best_k','recurrence_best_k','int'),('best_lag_bp','recurrence_best_lag_bp','int')])
    if len(hits)==0:
        row={'case':t['case'],'accession':t['accession'],'name':t['name'],'truth_start0':s0,'truth_end0':e0,'truth_type':t['truth_type'],
             'rediscovered':False,'clean_region_rank':'','clean_region_n':len(regs),'clean_region_rank_fraction':'','clean_region_id':'','clean_region_start0':'','clean_region_end0':'','support_family_n':'','support_families':'',
             'best_truth_acnv_global_rank':ba['global_rank'] if ba else '','best_truth_acnv_n_windows':ba['n_windows'] if ba else '',
             'best_truth_recurrence_global_rank':br['global_rank'] if br else '','best_truth_recurrence_n_windows':br['n_windows'] if br else '',
             'best_truth_recurrence_k':br.get('best_k','') if br else '','best_truth_recurrence_lag_bp':br.get('best_lag_bp','') if br else ''}
    else:
        # strict overlap merge means one truth interval can in principle overlap >1 clean regions; choose best rank but record count
        h=hits.sort_values('region_rank').iloc[0]
        row={'case':t['case'],'accession':t['accession'],'name':t['name'],'truth_start0':s0,'truth_end0':e0,'truth_type':t['truth_type'],
             'rediscovered':True,'truth_overlapping_clean_region_n':int(len(hits)),'clean_region_rank':int(h.region_rank),'clean_region_n':int(len(regs)),
             'clean_region_rank_fraction':float(h.region_rank/len(regs)),'clean_region_id':str(h.region_id),'clean_region_start0':int(h.start0),'clean_region_end0':int(h.end0),
             'support_family_n':int(h.support_family_n),'support_families':str(h.support_families),
             'best_q_acnv_family':float(h.best_q_acnv_family),'best_q_recurrence':float(h.best_q_recurrence),
             'recurrence_region_best_k':int(h.recurrence_best_k),'recurrence_region_best_lag_bp':int(h.recurrence_best_lag_bp),
             'best_truth_acnv_global_rank':ba['global_rank'] if ba else '','best_truth_acnv_n_windows':ba['n_windows'] if ba else '',
             'best_truth_recurrence_global_rank':br['global_rank'] if br else '','best_truth_recurrence_n_windows':br['n_windows'] if br else '',
             'best_truth_recurrence_k':br.get('best_k','') if br else '','best_truth_recurrence_lag_bp':br.get('best_lag_bp','') if br else ''}
    rows.append(row)
df=pd.DataFrame(rows);op=root/'outputs/CLEAN_TWO_ARM_POSTFREEZE_UNBLINDING.tsv';df.to_csv(op,sep='\t',index=False)
summary={'gate':'PASS_CLEAN_TWO_ARM_POSTFREEZE_REPLAY_COMPLETE','gate_pass':True,'scientific_status':'PAPER_SIMPLIFICATION_REPLAY_FROM_PREUNBLIND_FROZEN_WINDOW_SCORES',
         'n_cases':len(rows),'rediscovered_n':int(sum(bool(x['rediscovered']) for x in rows)),
         'clean_region_ranks':[x['clean_region_rank'] for x in rows],'clean_region_counts':[x['clean_region_n'] for x in rows],
         'cases':rows,'retuning_authorized':False,'sequence_features_recomputed':False,
         'interpretation':'Ranks are recomputed from the already-frozen pre-unblinding window scores using only ACNV_FAMILY and GENERIC_RECURRENCE with the predeclared q<=0.005 threshold, strict-overlap merging, and deterministic two-family ranking.'}
write_json(root/'outputs/CLEAN_TWO_ARM_FINAL_GATE.json',summary);print(json.dumps(summary,indent=2))
