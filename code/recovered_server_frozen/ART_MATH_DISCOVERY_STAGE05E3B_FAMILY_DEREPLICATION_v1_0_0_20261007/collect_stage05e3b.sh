#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
python - <<'PY'
import json
j=json.load(open('outputs/GATE_STAGE05E3B.json'))
assert j['gate_pass'] and j['gate']=='PASS_E2A_FAMILY_LINEAGE_DEREPLICATION_AND_CONTEXT_SYNTHESIS_FROZEN'
print('[PASS] Stage05E3B gate')
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="$HERE/ART_STAGE05E3B_FAMILY_DEREPLICATION_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" config scripts docs tests run_stage05e3b.sh collect_stage05e3b.sh outputs logs README.md
sha256sum "$OUT" > "$OUT.sha256"
echo "$OUT"
cat "$OUT.sha256"
