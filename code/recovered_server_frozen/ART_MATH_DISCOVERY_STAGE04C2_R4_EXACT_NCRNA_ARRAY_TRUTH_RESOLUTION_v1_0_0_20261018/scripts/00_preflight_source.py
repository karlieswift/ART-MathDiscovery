#!/usr/bin/env python3
from pathlib import Path
import argparse, json, sys
from common import sha256, md5, read_tsv, dump_json

EXPECTED_XLSX_MD5='75f92db2b03d1de60cb3958d28bd4021'
EXPECTED_R3_TRUTH_LOCK_SHA='ab694c29a1661b3c71cf2e600ed652b0f7e4bf0d66d7ae0b3bde0beb865aff19'

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ug-root',required=True)
    ap.add_argument('--r3-root',required=True)
    ap.add_argument('--out',required=True)
    a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    ug=Path(a.ug_root); r3=Path(a.r3_root)
    xlsx=ug/'inputs/Minerva_SupplementaryTables_preprint.xlsx'
    truth=r3/'outputs/UG27_R3_CASSETTE_TRUTH_LOCK.tsv'
    problems=[]
    if not xlsx.exists(): problems.append(f'missing workbook: {xlsx}')
    if not truth.exists(): problems.append(f'missing R3 cassette truth lock: {truth}')
    xmd5=md5(xlsx) if xlsx.exists() else None
    tsha=sha256(truth) if truth.exists() else None
    if xmd5 and xmd5 != EXPECTED_XLSX_MD5:
        problems.append(f'workbook MD5 mismatch: {xmd5}')
    if tsha and tsha != EXPECTED_R3_TRUTH_LOCK_SHA:
        problems.append(f'R3 truth lock SHA256 mismatch: {tsha}')
    rows=read_tsv(truth) if truth.exists() else []
    if rows and len(rows)!=5: problems.append(f'expected 5 truth rows, got {len(rows)}')
    gate={
      'stage':'Stage04C2-R4A-preflight',
      'gate':'PASS_STAGE04C2_R4A_AUTHOR_SOURCE_AND_R3_TRUTH_LOCK_PREFLIGHT' if not problems else 'HOLD_STAGE04C2_R4A_PREFLIGHT_FAILED',
      'gate_pass':not problems,
      'source_workbook_md5':xmd5,
      'expected_source_workbook_md5':EXPECTED_XLSX_MD5,
      'r3_cassette_truth_lock_sha256':tsha,
      'expected_r3_cassette_truth_lock_sha256':EXPECTED_R3_TRUTH_LOCK_SHA,
      'n_prespecified_systems':len(rows),
      'detector_outputs_read':False,
      'score_recomputed':False,
      'retuning_authorized':False,
      'science_change':False,
      'problems':problems
    }
    dump_json(out/'STAGE04C2_R4A_PREFLIGHT_GATE.json',gate)
    print(json.dumps(gate,indent=2))
    if problems: sys.exit(2)
if __name__=='__main__': main()
