#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
TS=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D4R4_PROSTT5_FOLDSEEK_REMOTE_STRUCTURE_RESULTS_${TS}.tar.gz"
FILES=(README.md config docs scripts tests run_selftest.sh prepare_foldseek_refs.sh run_stage05d4r4.sh collect_stage05d4r4.sh outputs logs raw/TIER_A_AFSP.m8)
tar -czf "$OUT" "${FILES[@]}"
sha256sum "$OUT" | tee "$OUT.sha256"
printf '%s\n' "$ROOT/$OUT"
