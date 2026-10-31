#!/usr/bin/env bash
set -euo pipefail
R5D=${1:?R5D_ROOT}
R5C=${2:?R5C_ROOT}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$R5D"
[ -f raw/P2_CLEAN_RT.m8 ] || { echo '[FAIL] existing raw/P2_CLEAN_RT.m8 missing'; exit 2; }
[ -s raw/P2_CLEAN_RT.m8 ] || { echo '[FAIL] existing raw/P2_CLEAN_RT.m8 empty'; exit 2; }
[ -f outputs/GATE_STAGE05D4R5D_PREDICTION.json ] || { echo '[FAIL] P2 prediction gate missing'; exit 2; }
python - <<'PY'
import json
x=json.load(open('outputs/GATE_STAGE05D4R5D_PREDICTION.json'))
assert x['gate_pass'] and x['n_P2']==40 and x['n_cif']==40 and x['n_confidence_json']==40,x
print('[PASS] existing P2 prediction gate 40/40')
PY
ROWS=$(wc -l < raw/P2_CLEAN_RT.m8)
echo "[PASS] reuse existing Foldseek m8 rows=$ROWS; Foldseek rerun=false"
python "$HERE/scripts/summarize_p2_v102.py" \
  --freeze frozen/STAGE05D4R5D_P2_FREEZE.tsv \
  --boltz boltz_p2_out \
  --m8 raw/P2_CLEAN_RT.m8 \
  --r5csummary "$R5C/outputs/STAGE05D4R5C_SUMMARY.json" \
  --outdir outputs
sha256sum \
  raw/P2_CLEAN_RT.m8 \
  outputs/STAGE05D4R5D_P2_CLEAN_RT_TRIAGE.tsv \
  outputs/STAGE05D4R5D_SUMMARY.json \
  outputs/GATE_STAGE05D4R5D.json \
  > outputs/STAGE05D4R5D_SHA256.tsv
cat outputs/STAGE05D4R5D_SUMMARY.json
echo '[PASS] Stage05D4-R5D v1.0.2 postprocessing closeout complete'
