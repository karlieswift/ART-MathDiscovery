#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd); stamp=$(date +%Y%m%d_%H%M%S); out="$S/ART_STAGE05F3V2_R1_60G_DELTA_RESULTS_${stamp}.tar.gz"
tar -czf "$out" -C "$S" frozen outputs logs --exclude='logs/per_accession/*' --exclude='outputs/formal_faa/*'
sha256sum "$out" | tee "$out.sha256"
echo "$out"
