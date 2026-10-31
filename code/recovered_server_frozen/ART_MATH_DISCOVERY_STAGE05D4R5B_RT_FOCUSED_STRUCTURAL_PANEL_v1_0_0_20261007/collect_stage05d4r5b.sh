#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
python - <<'PY'
import json
r=json.load(open('outputs/GATE_STAGE05D4R5B.json'))
assert r['gate_pass'] and r['gate']=='PASS_P1_RT_FOCUSED_STRUCTURAL_PANEL_SEARCH_COMPLETE',r
print('[PASS] R5B gate')
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D4R5B_RT_FOCUSED_STRUCTURAL_PANEL_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" config docs frozen outputs raw logs scripts README.md
sha256sum "$OUT" > "$OUT.sha256"
echo "[RESULT] $ROOT/$OUT"
cat "$OUT.sha256"
