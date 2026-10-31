#!/usr/bin/env python3
import argparse,json,tarfile
from pathlib import Path
import pandas as pd
from common import *
ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--marshill-root',required=True);ap.add_argument('--stage03a-tar',required=True);a=ap.parse_args()
root=Path(a.root).resolve(); mr=Path(a.marshill_root).resolve(); tar=Path(a.stage03a_tar).resolve(); pol=read_json(root/'config/TWO_ARM_POLICY.json');thr=float(pol['proposal_q_threshold'])
out=root/'outputs/boards';out.mkdir(parents=True,exist_ok=True)
manifest=[]
def write_case(case,d,source_hash,source_desc):
    cdir=out/case;cdir.mkdir(parents=True,exist_ok=True)
    sel,regs=build_two_arm_board(d,thr=thr,prefix='AR')
    sp=cdir/'CLEAN_TWO_ARM_SELECTED_PROPOSALS.tsv';rp=cdir/'CLEAN_TWO_ARM_UNION_REGIONS.tsv'
    sel.to_csv(sp,sep='\t',index=False);regs.to_csv(rp,sep='\t',index=False)
    gate={'gate':'PASS_CLEAN_TWO_ARM_BOARD_FROZEN','gate_pass':True,'case':case,'proposal_q_threshold':thr,
          'families':['ACNV_FAMILY','GENERIC_RECURRENCE'],'n_genomic_windows':int(len(d)),'n_selected_proposals':int(len(sel)),
          'n_union_regions':int(len(regs)),'source_frozen_scores_sha256':source_hash,'source':source_desc,
          'truth_used':False,'retuning_authorized':False,'recompute_sequence_features':False,
          'selected_sha256':sha256_file(sp),'union_regions_sha256':sha256_file(rp)}
    write_json(cdir/'CLEAN_TWO_ARM_BOARD_GATE.json',gate);manifest.append(gate);print(json.dumps(gate,indent=2))
ms=mr/'outputs/STAGE02B_WINDOW_FAMILY_SCORES.tsv.gz';d=pd.read_csv(ms,sep='\t');write_case('MW248466.1_MarsHill',d,sha256_file(ms),str(ms))
with tarfile.open(tar,'r:gz') as tf:
    for acc in ['MW349129','MW349128','MW218148.1']:
        d,member=read_tsv_from_tar(tf,f'outputs/{acc}/STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz')
        # hash exact compressed tar member bytes for provenance
        raw=tf.extractfile(member).read(); import hashlib; h=hashlib.sha256(raw).hexdigest()
        write_case(acc,d,h,f'{tar}::{member}')
write_json(root/'outputs/CLEAN_TWO_ARM_BOARD_MANIFEST.json',{'gate':'PASS_ALL_FOUR_CLEAN_TWO_ARM_BOARDS_FROZEN','gate_pass':True,'cases':manifest,'truth_used':False})
