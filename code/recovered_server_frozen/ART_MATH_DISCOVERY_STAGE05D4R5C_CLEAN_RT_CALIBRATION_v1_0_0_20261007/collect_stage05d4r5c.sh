#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd);cd "$ROOT"
python - <<'PY'
import json
g=json.load(open('outputs/GATE_STAGE05D4R5C.json'));assert g['gate_pass'] and g['gate']=='PASS_CLEAN_RT_STRUCTURAL_CALIBRATION_COMPLETE',g
print('[PASS] Stage05D4-R5C gate')
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D4R5C_CLEAN_RT_CALIBRATION_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" config docs frozen outputs raw logs scripts README.md
sha256sum "$OUT" > "$OUT.sha256"
echo "[RESULT] $ROOT/$OUT";cat "$OUT.sha256"
