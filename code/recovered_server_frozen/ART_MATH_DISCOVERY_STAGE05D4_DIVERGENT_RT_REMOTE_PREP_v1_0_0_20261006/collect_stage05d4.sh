#!/usr/bin/env bash
set -euo pipefail
TS=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D4_DIVERGENT_RT_REMOTE_PREP_RESULTS_${TS}.tar.gz"
tar -czf "$OUT" outputs logs config docs scripts tests README.md run_selftest.sh run_stage05d4.sh collect_stage05d4.sh
sha256sum "$OUT" | tee "$OUT.sha256"
echo "[DONE] $OUT"
