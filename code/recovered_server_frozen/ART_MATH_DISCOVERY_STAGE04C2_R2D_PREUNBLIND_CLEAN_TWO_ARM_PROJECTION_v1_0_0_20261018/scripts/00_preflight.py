#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, json, platform, sys
from pathlib import Path

EXPECTED = {
    'outputs/UG27_SCORE_FREEZE_GATE.json':'14f29cc2a856494308ffcce9b6c0d883eb10c988592ab8ebbac0e9d39447b498',
    'outputs/TWO_ARM_REGION_VIEW.tsv':'1a48f06cb2991f649502d843a0e853d011c51e45cf8ba3107a467e451d1e17f3',
    'outputs/HISTORICAL_EQUIVALENCE_GATE.json':'fe30d3519b693037c6e5d6705ccf06ef7129ffbac2379be09fe717518306c2ca',
    'outputs/SCORE_OUTPUT_MANIFEST.tsv':'fc5baa39e6108b87f21e8e96c4691e33d57741532fdbe01b3cd6d4955d9efc95',
}
HIST_EXPECTED = {
    'scripts/common.py':'a59b0e8439747487238b25249116f476baf54e55a5a6681411cb2272763cf255',
    'scripts/01_build_boards.py':'e0f7de502314f1854b2520b97892b56052b05b9db1f0be88c2ba15aa3493680a',
    'tests/test_two_arm.py':'8ff3f567e3eb6c040a105e5ca45038b9e840100365a3853430091812c7b902da',
    'README.md':'91ef9fef928a0969fa5393903874f5bbd3ee6bb81ac910add19246c5c8a0dee2',
}

def sha(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def fail(msg): raise SystemExit('HOLD_PREUNBLIND_PREFLIGHT: '+msg)

ap=argparse.ArgumentParser()
ap.add_argument('--root',required=True)
ap.add_argument('--r2c',required=True)
ap.add_argument('--historical-clean',required=True)
a=ap.parse_args()
root=Path(a.root).resolve(); r2c=Path(a.r2c).resolve(); hist=Path(a.historical_clean).resolve()
out=root/'outputs'; out.mkdir(parents=True,exist_ok=True)

for rel,h in EXPECTED.items():
    p=r2c/rel
    if not p.is_file(): fail(f'missing frozen R2C file {p}')
    got=sha(p)
    if got!=h: fail(f'R2C hash mismatch {rel}: {got} != {h}')

for rel,h in HIST_EXPECTED.items():
    p=hist/rel
    if not p.is_file(): fail(f'missing historical clean source {p}')
    got=sha(p)
    if got!=h: fail(f'historical implementation hash mismatch {rel}: {got} != {h}')

policy=hist/'config/TWO_ARM_POLICY.json'
if not policy.is_file(): fail(f'missing historical policy {policy}')
pol=json.loads(policy.read_text())
if float(pol.get('proposal_q_threshold',-1)) != 0.005:
    fail(f'historical proposal_q_threshold is not 0.005: {pol.get("proposal_q_threshold")}')

score_gate=json.loads((r2c/'outputs/UG27_SCORE_FREEZE_GATE.json').read_text())
checks={
 'gate_pass': score_gate.get('gate_pass') is True,
 'gate': score_gate.get('gate')=='PASS_ALL_15_CONTIG_BOARDS_AND_TWO_ARM_VIEW_FROZEN',
 'n_scoring_contigs': int(score_gate.get('n_scoring_contigs',-1))==15,
 'n_completed': int(score_gate.get('n_completed',-1))==15,
 'truth_used_for_scoring_false': score_gate.get('truth_used_for_scoring') is False,
 'truth_coordinates_loaded_false': score_gate.get('truth_coordinates_loaded') is False,
 'retuning_authorized_false': score_gate.get('retuning_authorized') is False,
}
if not all(checks.values()): fail('R2C score gate semantic checks failed: '+repr(checks))

# Parse frozen output manifest and verify all 15 exact window-score files.
man=r2c/'outputs/SCORE_OUTPUT_MANIFEST.tsv'
rows=list(csv.DictReader(man.open(),delimiter='\t'))
score_rows=[x for x in rows if x['path'].endswith('/STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz')]
if len(score_rows)!=15: fail(f'expected 15 frozen score files in manifest, got {len(score_rows)}')
for x in score_rows:
    p=r2c/'outputs'/x['path']
    if not p.is_file(): fail(f'missing frozen score file {p}')
    if p.stat().st_size!=int(x['bytes']): fail(f'byte-size mismatch {x["path"]}')
    got=sha(p)
    if got!=x['sha256']: fail(f'score hash mismatch {x["path"]}: {got} != {x["sha256"]}')

# Input lock must contain the same 15 scoring IDs.
lock=list(csv.DictReader((r2c/'outputs/CONTIG_INPUT_LOCK.tsv').open(),delimiter='\t'))
if len(lock)!=15: fail(f'expected 15 CONTIG_INPUT_LOCK rows, got {len(lock)}')
ids={x['scoring_id'] for x in lock}
manifest_ids={Path(x['path']).parts[1] for x in score_rows}
if ids!=manifest_ids: fail('CONTIG_INPUT_LOCK scoring IDs != score manifest scoring IDs')

# Guard: this package never requires/opens R3 truth.
truth_like=[]
for p in root.rglob('*'):
    if p.is_file() and any(t in p.name.upper() for t in ['TRUTH','CONSTRUCT_NATIVE_BRIDGE','POSTFREEZE_UNBLIND']):
        truth_like.append(str(p.relative_to(root)))
if truth_like: fail('unexpected truth-like files inside pre-unblind package: '+repr(truth_like))

pre={
 'stage':'Stage04C2-R2D-preflight',
 'gate':'PASS_STAGE04C2_R2D_PREUNBLIND_INPUTS_AND_HISTORICAL_CLEAN_RULE_LOCKED',
 'gate_pass':True,
 'r2c_root':str(r2c),
 'historical_clean_root':str(hist),
 'n_scoring_contigs':15,
 'frozen_score_files_verified_n':15,
 'historical_proposal_q_threshold':0.005,
 'historical_policy_sha256':sha(policy),
 'historical_implementation_hashes':HIST_EXPECTED,
 'r2c_principal_hashes':EXPECTED,
 'truth_loaded':False,
 'truth_paths_read':False,
 'retuning_authorized':False,
 'score_recomputed':False,
 'q_recomputed':False,
 'sequence_features_recomputed':False,
 'science_change':False,
 'runtime':{'python':sys.version,'platform':platform.platform()},
}
(out/'STAGE04C2_R2D_PREFLIGHT_GATE.json').write_text(json.dumps(pre,indent=2,sort_keys=True)+'\n')
print(json.dumps(pre,indent=2,sort_keys=True))
