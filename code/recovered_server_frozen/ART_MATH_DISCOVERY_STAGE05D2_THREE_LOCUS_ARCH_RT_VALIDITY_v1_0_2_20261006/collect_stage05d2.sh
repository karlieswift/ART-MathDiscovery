#!/usr/bin/env bash
set -euo pipefail
TS=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D2_THREE_LOCUS_RESULTS_${TS}.tar.gz"
tar -czf "$OUT" README.md config docs scripts outputs work/ncbi 2>/dev/null || tar -czf "$OUT" README.md config docs scripts outputs
sha256sum "$OUT" > "$OUT.sha256"
echo "$OUT"
echo "$OUT.sha256"
