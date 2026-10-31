#!/usr/bin/env bash
set -euo pipefail
TS=$(date +%Y%m%d_%H%M%S)
ARCH="ART_STAGE05D4R3_COMPACT_STRUCTURE_READY_RESULTS_${TS}.tar.gz"
tar -czf "$ARCH" outputs logs config docs scripts tests README.md run_selftest.sh run_stage05d4r3.sh collect_stage05d4r3.sh
sha256sum "$ARCH" | tee "$ARCH.sha256"
echo "$PWD/$ARCH"
