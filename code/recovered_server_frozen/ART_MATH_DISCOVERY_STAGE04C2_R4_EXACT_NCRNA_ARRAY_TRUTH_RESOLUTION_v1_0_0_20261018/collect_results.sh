#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_MATH_DISCOVERY_STAGE04C2_R4_EXACT_NCRNA_ARRAY_TRUTH_RESULTS_${TS}.tar.gz"
tar -czf "$PKG" outputs logs README.md config docs PACKAGE_MANIFEST.sha256 run_selftest.sh run_stage04c2_r4a_truth_resolution.sh run_stage04c2_r4b_exact_array_eval.sh scripts tests
sha256sum "$PKG" > "${PKG}.sha256"
sha256sum -c "${PKG}.sha256"
ls -lh "$PKG" "${PKG}.sha256"
