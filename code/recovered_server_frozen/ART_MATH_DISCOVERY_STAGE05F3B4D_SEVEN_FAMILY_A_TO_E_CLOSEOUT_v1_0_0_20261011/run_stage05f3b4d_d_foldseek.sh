#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
RTDB="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007/references/clean_rt_afsp"
BSTAT="$S/outputs/F3B4D_D_BOLTZ_STATUS.tsv"
[[ -s "$BSTAT" ]] || { echo '[FATAL] run D Boltz first'; exit 2; }
command -v foldseek >/dev/null || { echo '[FATAL] foldseek unavailable; activate art_novelty'; exit 2; }
ls "${RTDB}"* >/dev/null 2>&1 || { echo '[FATAL] frozen R5C DB missing'; exit 2; }
python "$S/scripts/run_d_foldseek.py" "$BSTAT" "$RTDB" "$S/raw/foldseek_d" "$S/tmp/foldseek_d" "$S/logs" "$THREADS" | tee "$S/logs/d_foldseek_summary.log"
python "$S/scripts/summarize_d.py" \
  "$S/outputs/F3B4D_D_STRUCTURE_EXTENSION_TARGETS_6.tsv" "$BSTAT" "$S/outputs/F3B4D_D_FOLDSEEK_STATUS.tsv" "$S/raw/foldseek_d" "$S/outputs" "$S/frozen/F3B4D_IMPLEMENTATION_LOCK.json" \
  2>&1 | tee "$S/logs/d_structure_summary.log"
echo PASS_STAGE05F3B4D_D_FOLDSEEK_AND_R5C_COMPLETE
