#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
test -f outputs/STAGE04B5_CLOSEOUT_GATE.json || { echo 'HOLD no final closeout gate'; exit 2; }
python - <<'PY'
import json
x=json.load(open('outputs/STAGE04B5_CLOSEOUT_GATE.json'))
assert x['gate_pass'] is True
assert x['gate']=='PASS_STAGE04B5_54G_MATCHED_ARCHITECTURE_REPLAY_COMPLETE'
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B5_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" \
  README.md config docs logs \
  outputs/STAGE04B5_CLOSEOUT_GATE.json \
  outputs/CONTEXT_EVALUABILITY_SUMMARY.json \
  outputs/ARCHITECTURE_BY_ARM.tsv \
  outputs/ARCHITECTURE_SELECTED_FAMILIES.tsv \
  outputs/DELTA26396_FAMILY_ARCHITECTURE_SUMMARY.tsv \
  outputs/DELTA88487_CANDIDATE_FAMILY_MAP.tsv \
  outputs/DELTA88487_CONTEXT_EVALUABILITY.tsv \
  outputs/DELTA88487_REPEAT_ARRAYS.tsv \
  outputs/DELTA88487_NEAREST_DOWNSTREAM_ORFS_WITH_CLUSTER.tsv \
  outputs/SHA256SUMS.txt
sha256sum "$OUT" | tee "$OUT.sha256"
echo "RESULT=$OUT"
