#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, shutil, tarfile, tempfile
from pathlib import Path

def sha256(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--archive',required=True)
    a=ap.parse_args()
    here=Path(__file__).resolve().parents[1]
    arc=Path(a.archive).resolve()
    if not arc.is_file(): raise SystemExit(f'archive not found: {arc}')
    arcsha=sha256(arc)
    with tempfile.TemporaryDirectory(prefix='f4a_b1_import_') as td:
        td=Path(td)
        with tarfile.open(arc,'r:gz') as tf: tf.extractall(td)
        gatep=td/'outputs/score_freeze/BASELINE_B1_GATE.json'
        seqp=td/'frozen/BASELINE_SEQUENCE_MANIFEST.tsv'
        invp=td/'frozen/BASELINE_BOARD_INVENTORY.tsv'
        if not (gatep.is_file() and seqp.is_file() and invp.is_file()):
            raise SystemExit('archive lacks required R2 B1 frozen state')
        gate=json.loads(gatep.read_text())
        if not gate.get('gate_pass'): raise SystemExit('archived B1 gate is not PASS')
        if gate.get('truth_opened_for_scoring') is not False: raise SystemExit('archived B1 chronology invalid')
        if 'R2' not in str(gate.get('gate','')): raise SystemExit('expected R2 B1 archive')
        # Replace only imported state; R4 config/scripts remain untouched.
        for rel in ['frozen','outputs/score_freeze']:
            src=td/rel; dst=here/rel
            if dst.exists(): shutil.rmtree(dst)
            shutil.copytree(src,dst,copy_function=shutil.copy2)
    obj={
      'stage':'Stage05F4A-B2-R4-import-frozen-B1',
      'version':'v1.0.4',
      'gate':'PASS_STAGE05F4A_B2_R4_R2_B1_FROZEN_STATE_IMPORTED',
      'gate_pass':True,
      'source_archive':str(arc),
      'source_archive_sha256':arcsha,
      'source_b1_gate':gate,
      'truth_opened_for_scoring':False,
      'retuning_authorized':False,
      'note':'Imported frozen B0/B1 artifacts byte-for-byte from the archived R2 score-freeze; B0/B1 are not rerun in R4.'
    }
    (here/'frozen').mkdir(exist_ok=True)
    (here/'frozen/IMPORTED_R2_B1_ARCHIVE.json').write_text(json.dumps(obj,indent=2)+"\n")
    print(json.dumps(obj,indent=2))
if __name__=='__main__': main()
