#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="$S/ART_STAGE05F3B4E_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" -C "$S" README_CN.txt frozen inputs outputs logs
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$OUT"
