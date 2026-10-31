#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path
from common import *
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); a=ap.parse_args(); root=Path(a.root).resolve()
    pol=read_json(root/'config/ENGINE_POLICY.json')
    exp={}
    with open(root/'config/FROZEN_SOURCE_SHA256.tsv') as f:
        for r in csv.DictReader(f,delimiter='\t'): exp[r['filename']]=r['sha256']
    rec=[]; ok=True
    for fn,h in exp.items():
        p=root/'frozen_sources'/fn; ex=p.exists(); obs=sha256_file(p) if ex else None; match=(obs==h); ok &= bool(match); rec.append({'filename':fn,'expected_sha256':h,'observed_sha256':obs,'hash_ok':match})
    cohort=[]
    with open(root/'config/VALIDATION_COHORT.tsv') as f: cohort=list(csv.DictReader(f,delimiter='\t'))
    ok &= len(cohort)==3 and all(x['stage02_development_used']=='NO' for x in cohort)
    gate={'gate':'PASS_STAGE03A_FROZEN_ENGINE_BEFORE_HOLDOUT_SCORING' if ok else 'FAIL_STAGE03A_FROZEN_ENGINE_BEFORE_HOLDOUT_SCORING','gate_pass':ok,'scientific_status':pol['scientific_status'],'blindness_level':pol['blindness_level'],'truth_policy':pol['truth_policy'],'validation_accessions':[x['accession'] for x in cohort],'source_hashes':rec,'development_freeze':pol['development_freeze'],'authorization':'If PASS, validation FASTA may be fetched and scored. No ART truth/family labels may be queried until Stage03A blind outputs are frozen.'}
    write_json(root/'outputs/STAGE03A_ENGINE_FREEZE_GATE.json',gate); print(json.dumps(gate,indent=2))
    if not ok: raise SystemExit(2)
if __name__=='__main__': main()
