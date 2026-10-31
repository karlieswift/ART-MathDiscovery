#!/usr/bin/env bash
set -euo pipefail
TS=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D3_REMAINING_BROAD_RT_RESULTS_${TS}.tar.gz"
tar -czf "$OUT" README.md config docs scripts outputs inputs/official_pc run_selftest.sh run_stage05d3.sh collect_stage05d3.sh windows_download_official_genbank.ps1
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$OUT"
