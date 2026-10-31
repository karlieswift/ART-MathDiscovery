#!/usr/bin/env python3
import csv, json, subprocess, sys, tempfile
from pathlib import Path

def main():
    root=Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory() as td:
        td=Path(td); src=td/'regions.tsv'; uni=td/'universe.tsv'; out=td/'out'
        cols=['analysis_id','sequence_id','arm','length_bp','region_rank','region_id','start0','end0','span_bp','selected_window_rows','best_q_acnv_family','best_q_recurrence','sum_family_surprisal','best_single_family_q','acnv_best_member','acnv_best_window_start0','acnv_best_window_end0','acnv_best_window_len','recurrence_best_k','recurrence_best_lag_bp','recurrence_best_window_start0','recurrence_best_window_end0','recurrence_best_window_len','support_families']
        rows=[
          ['A','A','L4_CURATED',1000,1,'R1',10,20,10,2,0.001,0.002,6,0.001,'X',0,100,100,8,50,0,100,100,'ACNV_FAMILY;GENERIC_RECURRENCE'],
          ['A','A','L4_CURATED',1000,2,'R2',30,40,10,1,0.001,0.2,3,0.001,'X',0,100,100,8,50,0,100,100,'ACNV_FAMILY'],
          ['B','B','L4_CURATED',2000,1,'R1',50,60,10,1,0.2,0.001,3,0.001,'X',0,100,100,8,50,0,100,100,'GENERIC_RECURRENCE']]
        with src.open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t'); w.writerow(cols); w.writerows(rows)
        with uni.open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t'); w.writerow(['analysis_id','sequence_id','arm','source','release_date','cluster_id','notes','length_bp','fasta_sha256','fasta_header','frozen_fasta']); w.writerow(['A','A','L4_CURATED','','','','',1000,'','','']); w.writerow(['B','B','L4_CURATED','','','','',2000,'','',''])
        subprocess.run([sys.executable,str(root/'scripts/build_stage04a_r1.py'),'--source',str(src),'--universe',str(uni),'--outdir',str(out),'--allow-noncanonical-source'],check=True,stdout=subprocess.DEVNULL)
        g=json.loads((out/'STAGE04A_R1_SUMMARY.json').read_text())
        assert g['n_all_frozen_union_regions']==3
        assert g['n_dual_arm_regions']==1
        assert g['n_sequences_with_ge1_dual_arm']==1
        assert g['score_recomputed'] is False
        print('[PASS] Stage04A-R1 two-arm projection selftest')
if __name__=='__main__': main()
