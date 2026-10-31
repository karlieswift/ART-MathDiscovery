#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
cd "$S"
[[ -s outputs/STAGE05F3B4D_GATE.json ]] || { echo '[FATAL] finalize first'; exit 2; }
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F3B4D_SEVEN_FAMILY_A_TO_E_CLOSEOUT_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" README_CN.txt frozen outputs logs inputs/F3B4D_D_STRUCTURE_EXTENSION_6.faa inputs/F3B4D_RT_ARCHITECTURE_REPRESENTATIVES_7.faa inputs/F3B4D_PARTNER_REPRESENTATIVES_7.faa structure_boltz_out raw/foldseek_d
sha256sum "$OUT" > "$OUT.sha256"
echo "$S/$OUT"
cat "$OUT.sha256"
