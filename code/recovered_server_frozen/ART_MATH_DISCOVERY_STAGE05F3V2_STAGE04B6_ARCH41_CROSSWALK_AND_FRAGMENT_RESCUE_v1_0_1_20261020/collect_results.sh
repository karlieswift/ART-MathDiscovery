#!/usr/bin/env bash
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd)
test -f "$H/outputs/STAGE04B6_CLOSEOUT_GATE.json" || { echo HOLD_NO_FINAL_GATE; exit 2; }
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="$H/ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B6_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" -C "$H" README.md config docs frozen outputs logs
sha256sum "$OUT" | tee "$OUT.sha256"
echo "RESULT=$OUT"
