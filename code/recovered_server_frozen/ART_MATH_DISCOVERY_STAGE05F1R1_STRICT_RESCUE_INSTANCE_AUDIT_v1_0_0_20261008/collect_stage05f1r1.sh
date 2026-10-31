#!/usr/bin/env bash
set -euo pipefail
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F1R1_STRICT_RESCUE_INSTANCE_AUDIT_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" outputs logs README.md docs
sha256sum "$OUT" > "$OUT.sha256"
echo "$PWD/$OUT"
cat "$OUT.sha256"
