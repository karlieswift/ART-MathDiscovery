#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
[ -f outputs/GATE_STAGE05E2.json ] || { echo '[FAIL] no E2 gate' >&2; exit 2; }
python - <<'PY'
import json
x=json.load(open('outputs/GATE_STAGE05E2.json'))
assert x.get('gate_pass') is True
assert x.get('n_E1_primary_families')==66
print('[PASS] E2 gate ready for collection')
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="../ART_STAGE05E2_GLOBAL_SPECIFICITY_TOPOLOGY_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" outputs config docs scripts README.md run_stage05e2.sh run_selftest.sh collect_stage05e2.sh
sha256sum "$OUT" > "$OUT.sha256"
echo "$OUT"
echo "$OUT.sha256"
