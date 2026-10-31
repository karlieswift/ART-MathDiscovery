#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
[[ -f "$HERE/outputs/STAGE04B8_FINAL_CLOSEOUT_GATE.json" ]] || { echo 'HOLD no final gate'; exit 2; }
python - "$HERE/outputs/STAGE04B8_FINAL_CLOSEOUT_GATE.json" <<'PY'
import json,sys
z=json.load(open(sys.argv[1])); assert z.get('gate_pass') and z.get('gate')=='PASS_STAGE04B8_FINAL_54G_CLOSEOUT'
PY
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="$HERE/../ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B8_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" -C "$HERE" outputs logs config docs README.md run_selftest.sh run_stage04b8_closeout.sh collect_results.sh scripts/stage04b8.py
sha256sum "$OUT" > "$OUT.sha256"
echo "$OUT"
echo "$OUT.sha256"
