#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
cd "$S"
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_STAGE05F3B4A_SURVIVOR_NATIVE_LOCUS_RECONSTRUCTION_RESULTS_${TS}.tar.gz"
tar -czf "$PKG" outputs logs frozen README_CN.txt run_stage05f3b4a.sh collect_stage05f3b4a.sh scripts tests
sha256sum "$PKG" | tee "$PKG.sha256"
ls -lh "$PKG" "$PKG.sha256"
