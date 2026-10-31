#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_STAGE05F3B4B1_FAMILY_WIDE_REPEAT_ARRAY_ARCHITECTURE_RESULTS_${TS}.tar.gz"
cd "$HERE"
tar -czf "$PKG" frozen outputs logs scripts tests README_CN.txt run_stage05f3b4b1.sh collect_stage05f3b4b1.sh
sha256sum "$PKG" > "$PKG.sha256"
echo "$HERE/$PKG"
cat "$PKG.sha256"
