#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, io, json, tarfile
from pathlib import Path
from common import sha256_file

EXPECTED_ARCHIVE_SHA='92fb006869bd868f9ced428ec23ad0f8e323cf04bc76d5454093a7b856dba46c'
EXPECTED_R2D={
 'outputs/STAGE04C2_R2D_GATE.json':'ecfa7e3b6a6418db660182fdbbd6e82135a39083be10ab4779c3ae581487cb6b',
 'outputs/CLEAN_TWO_ARM_ALL_REGIONS.tsv':'1d8f394a5e35e48301c4c0ebe959b2babb301269eb5977bb9c66c9c86545e3c9',
 'outputs/CLEAN_TWO_ARM_DUAL_ARM_VIEW.tsv':'12940178a8efbff191c066e9bec28ed60ce1b4dad1da7507ba9049995c9dbe15',
 'outputs/CLEAN_VS_LEGACY_DUAL_ARM_SCORE_ONLY_COMPARISON.tsv':'7bc7b212e6e889a007a9a474f39049a708f5a788b6d3c92230145e4c4be411d0',
}
EXPECTED_R2C_INPUT_LOCK='2d60695a3ab7ce0b0b09ff3a68236201789918f78ff29c035f3bbe56fa867f21'

def fail(x): raise SystemExit('HOLD_STAGE04C2_R3_SCORE_ARCHIVE_PREFLIGHT: '+x)

def sha_bytes(b): return hashlib.sha256(b).hexdigest()

ap=argparse.ArgumentParser()
ap.add_argument('--root',required=True)
ap.add_argument('--r2d',required=True)
ap.add_argument('--r2c',required=True)
ap.add_argument('--r2d-archive',required=True)
a=ap.parse_args()
root=Path(a.root).resolve(); r2d=Path(a.r2d).resolve(); r2c=Path(a.r2c).resolve(); arc=Path(a.r2d_archive).resolve()
out=root/'outputs'; out.mkdir(parents=True,exist_ok=True)

if not arc.is_file(): fail(f'missing R2D frozen result archive {arc}')
if sha256_file(arc)!=EXPECTED_ARCHIVE_SHA: fail('R2D result archive SHA256 mismatch')

for rel,h in EXPECTED_R2D.items():
    p=r2d/rel
    if not p.is_file(): fail(f'missing live R2D output {p}')
    if sha256_file(p)!=h: fail(f'live R2D hash mismatch {rel}')

inp=r2c/'outputs/CONTIG_INPUT_LOCK.tsv'
if not inp.is_file() or sha256_file(inp)!=EXPECTED_R2C_INPUT_LOCK:
    fail('R2C CONTIG_INPUT_LOCK hash mismatch')

# Verify the archived R2D tar itself contains the same exact frozen boards/gate.
with tarfile.open(arc,'r:gz') as tf:
    names=tf.getnames()
    for rel,h in EXPECTED_R2D.items():
        hits=[n for n in names if n.endswith(rel)]
        if len(hits)!=1: fail(f'archive member ambiguity for {rel}: {hits}')
        raw=tf.extractfile(hits[0]).read()
        if sha_bytes(raw)!=h: fail(f'archived R2D member hash mismatch {rel}')

G=json.loads((r2d/'outputs/STAGE04C2_R2D_GATE.json').read_text())
checks={
 'gate':G.get('gate')=='PASS_STAGE04C2_R2D_PREUNBLIND_CLEAN_ACNV_RECURRENCE_BOARD_FROZEN',
 'gate_pass':G.get('gate_pass') is True,
 'n_scoring_contigs':int(G.get('n_scoring_contigs',-1))==15,
 'n_clean_union_regions':int(G.get('n_clean_union_regions',-1))==6212,
 'n_clean_dual_arm_regions':int(G.get('n_clean_dual_arm_regions',-1))==303,
 'truth_loaded_false':G.get('truth_loaded') is False,
 'truth_paths_read_false':G.get('truth_paths_read') is False,
 'cassette_overlap_false':G.get('cassette_overlap_computed') is False,
 'retuning_false':G.get('retuning_authorized') is False,
 'score_recomputed_false':G.get('score_recomputed') is False,
 'q_recomputed_false':G.get('q_recomputed') is False,
 'sequence_features_recomputed_false':G.get('sequence_features_recomputed') is False,
 'science_change_false':G.get('science_change') is False,
}
if not all(checks.values()): fail('R2D semantic gate checks failed: '+repr(checks))

D={
 'stage':'Stage04C2-R3-score-archive-preflight',
 'gate':'PASS_STAGE04C2_R3_R2D_ARCHIVE_AND_LIVE_CLEAN_BOARD_BYTE_LOCKED_BEFORE_TRUTH_READ',
 'gate_pass':True,
 'r2d_result_archive':str(arc),
 'r2d_result_archive_sha256':EXPECTED_ARCHIVE_SHA,
 'r2d_live_hashes':EXPECTED_R2D,
 'r2c_contig_input_lock_sha256':EXPECTED_R2C_INPUT_LOCK,
 'truth_opened':False,
 'cassette_coordinates_read':False,
 'retuning_authorized':False,
 'score_recomputed':False,
 'q_recomputed':False,
 'science_change':False,
}
(out/'STAGE04C2_R3_SCORE_ARCHIVE_PREFLIGHT_GATE.json').write_text(json.dumps(D,indent=2,sort_keys=True)+'\n')
print(json.dumps(D,indent=2,sort_keys=True))
