#!/usr/bin/env python3
from pathlib import Path
import argparse, json, subprocess, sys, zipfile
from common import sha256_file, write_json

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); a=ap.parse_args()
    root=Path(a.root).resolve(); pol=json.loads((root/'config/STAGE04A_POLICY.json').read_text())
    z=root/'provenance'/pol['stage03a_zip']
    got=sha256_file(z)
    if got != pol['stage03a_zip_sha256']:
        raise SystemExit(f'Stage03A ZIP SHA256 mismatch: {got}')
    # Stage04A must not contain Stage03B truth files.
    prohibited_names={'ART_TRUTH_MANIFEST.tsv','ART_REPEAT_ANCHORS.tsv'}
    found=[]
    for p in root.rglob('*'):
        if p.is_file() and p.name in prohibited_names: found.append(str(p.relative_to(root)))
    if found: raise SystemExit('Prohibited truth files present in Stage04A package: '+','.join(found))
    runtime=root/'runtime_engine'
    if runtime.exists():
        import shutil; shutil.rmtree(runtime)
    runtime.mkdir(parents=True)
    with zipfile.ZipFile(z) as zz: zz.extractall(runtime)
    kids=[p for p in runtime.iterdir() if p.is_dir()]
    if len(kids)!=1: raise SystemExit('unexpected Stage03A zip layout')
    eng=kids[0]
    # Byte parity of frozen scientific sources is checked by Stage03A selftest.
    cp=subprocess.run(['bash',str(eng/'run_selftest.sh')],cwd=eng,text=True,capture_output=True)
    (root/'logs').mkdir(exist_ok=True)
    (root/'logs'/'00_stage03a_selftest.log').write_text(cp.stdout+'\n'+cp.stderr)
    if cp.returncode!=0: raise SystemExit('embedded Stage03A selftest failed; see logs/00_stage03a_selftest.log')
    gate={
      'gate':'PASS_STAGE04A_BYTE_EXACT_STAGE03A_ENGINE_PREFLIGHT',
      'gate_pass':True,
      'stage03a_zip_sha256':got,
      'runtime_engine_root':str(eng.relative_to(root)),
      'truth_files_present':False,
      'retuning_authorized':False
    }
    write_json(root/'outputs/STAGE04A_PREFLIGHT_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
