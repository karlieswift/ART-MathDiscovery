#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B4R1_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" README.md config docs scripts run_stage04b4r1.sh collect_results.sh outputs logs PACKAGE_SHA256.tsv
sha256sum "$OUT" | tee "$OUT.sha256"
echo "RESULT=$OUT"
