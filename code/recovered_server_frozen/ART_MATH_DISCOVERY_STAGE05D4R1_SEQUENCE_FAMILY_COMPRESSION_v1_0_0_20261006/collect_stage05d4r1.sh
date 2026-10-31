#!/usr/bin/env bash
set -euo pipefail
TS=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D4R1_SEQUENCE_FAMILY_COMPRESSION_RESULTS_${TS}.tar.gz"
tar -czf "$OUT" outputs logs config docs scripts tests README.md run_selftest.sh run_stage05d4r1.sh collect_stage05d4r1.sh
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$OUT"
