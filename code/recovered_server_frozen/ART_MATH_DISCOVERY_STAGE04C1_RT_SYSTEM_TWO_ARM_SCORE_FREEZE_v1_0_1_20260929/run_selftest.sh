#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
python - "$ROOT" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]); pol=json.loads((root/'config/C1_POLICY.json').read_text())
assert pol['version']=='v1.0.1'
assert pol['eligible_control_ids']==['SEN2_14028S','EC78_ECO7','DRT2_KPN']
assert pol['proposal_q_threshold']==0.005
assert pol['recurrence']=={'k_values':[4,6,8],'min_lag_bp':8,'max_lag_bp':512,'min_valid_pairs':32,'min_recurrence_match_count':2,'truth_specific_lag_band':False}
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(root/'frozen_sources/ACNV_SOURCE__legacy_acnv_extractor_frozen.py')=='2e5754476d7615044dcb50d131e45e3f18b214aea5540d39f5bcb05ceb2d2a28'
assert sha(root/'scripts/common_stage03a_exact.py')=='d215ddebdf776bfa4ea7e0ef362c2ccb6407278df4bf1d4ebbc87dfd68f95424'
txt=(root/'run_stage04c1_score_freeze.sh').read_text()+(root/'scripts/00_preflight.py').read_text()
assert 'STAGE04C0_R3_FIXED_SEQUENCE_MANIFEST.tsv' in txt
assert 'STAGE04C0_R3_FINAL_GATE.json' in txt
assert 'STAGE04C0_R1_FIXED_SEQUENCE_MANIFEST.tsv' not in txt
print('[PASS] Stage04C1 v1.0.1 frozen two-arm policy/source/R3-input selftest')
PY
