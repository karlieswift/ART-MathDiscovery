#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
python - <<'PY'
import json
for p in ['frozen/STAGE05D4R5D_FREEZE.json','outputs/GATE_STAGE05D4R5D_PREDICTION.json','outputs/GATE_STAGE05D4R5D.json']:
 x=json.load(open(p));assert x['gate_pass'],(p,x)
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D4R5D_P2_CALIBRATED_FULL3D_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" config docs frozen outputs raw logs scripts README.md run_selftest.sh prepare_stage05d4r5d.sh run_stage05d4r5d_p2.sh run_clean_rt_screen.sh collect_stage05d4r5d.sh
sha256sum "$OUT" > "$OUT.sha256"
echo "[PASS] $OUT"
cat "$OUT.sha256"
