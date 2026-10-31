#!/usr/bin/env python3
import argparse,json
from pathlib import Path
from common import *
ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--marshill-root',required=True);ap.add_argument('--stage03a-tar',required=True);a=ap.parse_args()
root=Path(a.root).resolve(); mr=Path(a.marshill_root).resolve(); tar=Path(a.stage03a_tar).resolve(); pol=read_json(root/'config/TWO_ARM_POLICY.json')
ms=mr/'outputs/STAGE02B_WINDOW_FAMILY_SCORES.tsv.gz'
if not ms.is_file(): raise SystemExit(f'FAIL missing MarsHill frozen score file: {ms}')
if not tar.is_file(): raise SystemExit(f'FAIL missing Stage03A results tar: {tar}')
th=sha256_file(tar); exp=pol['expected_stage03a_results_tar_sha256']
if th!=exp: raise SystemExit(f'FAIL Stage03A tar SHA256 mismatch\nexpected={exp}\nactual  ={th}')
# optional historical gate cross-check if it exposes the score hash
mg=mr/'outputs/STAGE02B_UNION_BOARD_GATE.json'; verified=None
if mg.is_file():
    try:
        g=read_json(mg)
        actual=sha256_file(ms)
        for k in ['window_family_scores_sha256','family_scores_sha256','window_scores_sha256']:
            if k in g:
                if str(g[k])!=actual: raise SystemExit(f'FAIL MarsHill score hash mismatch vs {k}')
                verified=k
                break
    except json.JSONDecodeError: pass
obj={'gate':'PASS_CLEAN_TWO_ARM_PREFLIGHT_FROZEN_INPUTS','gate_pass':True,
     'scientific_status':pol['scientific_status'],'stage03a_results_tar':str(tar),'stage03a_results_tar_sha256':th,
     'marshill_score_file':str(ms),'marshill_score_sha256':sha256_file(ms),'marshill_historical_gate_hash_check_key':verified,
     'recompute_sequence_features':False,'truth_used':False,'retuning_authorized':False}
write_json(root/'outputs/CLEAN_TWO_ARM_PREFLIGHT_GATE.json',obj);print(json.dumps(obj,indent=2))
