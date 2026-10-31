#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
TS=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F4A_FORMAL_BASELINE_RESULTS_${TS}.tar.gz"
tar -czf "$OUT" config frozen outputs logs README_CN.md 00_preflight_and_freeze.sh 01_run_score_freeze.sh 02_postfreeze_unblind.sh scripts tests
sha256sum "$OUT" > "${OUT}.sha256"
echo "$HERE/$OUT"
cat "${OUT}.sha256"
