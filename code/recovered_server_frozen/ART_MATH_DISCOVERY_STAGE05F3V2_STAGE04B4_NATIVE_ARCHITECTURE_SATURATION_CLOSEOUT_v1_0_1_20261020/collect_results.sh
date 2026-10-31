#!/usr/bin/env bash
set -euo pipefail
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B4_RESULTS_${TS}.tar.gz"
tar -czf "$PKG" outputs logs frozen config docs scripts README.md run_selftest.sh run_stage04b4_closeout.sh collect_results.sh PACKAGE_SHA256.tsv
sha256sum "$PKG" | tee "${PKG}.sha256"
echo "RESULT=$PKG"
