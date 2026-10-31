#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"; cd "$ROOT"
TS=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D4R2_SWISSPROT_SENSITIVE_HOMOLOGY_RESULTS_${TS}.tar.gz"
tar -czf "$OUT" outputs logs config docs scripts tests README.md run_selftest.sh run_stage05d4r2.sh collect_stage05d4r2.sh
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$ROOT/$OUT"
