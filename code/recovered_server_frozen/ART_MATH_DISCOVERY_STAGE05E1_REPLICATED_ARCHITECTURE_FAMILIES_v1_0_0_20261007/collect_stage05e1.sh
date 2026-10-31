#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
python - <<'PY'
import json
j=json.load(open('outputs/GATE_STAGE05E1.json'))
assert j['gate_pass'] and j['gate']=='PASS_REPLICATED_ARCHITECTURE_FAMILIES_FROZEN'
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05E1_REPLICATED_ARCHITECTURE_FAMILIES_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" README.md config docs scripts run_selftest.sh run_stage05e1.sh collect_stage05e1.sh outputs logs
sha256sum "$OUT" > "$OUT.sha256"
echo "[PASS] $OUT"
cat "$OUT.sha256"
