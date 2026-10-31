#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F3B4B0_FAMILY_WIDE_CONTEXT_COMPLETENESS_RESCUE_RESULTS_${STAMP}.tar.gz"
cd "$HERE"
tar -czf "$OUT" outputs frozen logs README_CN.txt scripts run_stage05f3b4b0.sh tests
sha256sum "$OUT" > "$OUT.sha256"
echo "$HERE/$OUT"
echo "$HERE/$OUT.sha256"
