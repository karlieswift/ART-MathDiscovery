#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
TS=$(date +%Y%m%d_%H%M%S)
OUT="$HERE/ART_STAGE05F3B4B1D_TECHNICAL_RESCUE_REPEAT_EXTENSION_RESULTS_${TS}.tar.gz"
tar -czf "$OUT" -C "$HERE" frozen outputs logs README_CN.txt run_stage05f3b4b1d.sh scripts tests
sha256sum "$OUT" > "$OUT.sha256"
ls -lh "$OUT" "$OUT.sha256"
cat "$OUT.sha256"
