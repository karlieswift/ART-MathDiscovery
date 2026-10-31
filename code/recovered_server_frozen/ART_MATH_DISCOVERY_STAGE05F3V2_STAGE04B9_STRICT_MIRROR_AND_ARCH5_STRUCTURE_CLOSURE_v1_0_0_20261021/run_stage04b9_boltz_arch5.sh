#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root}; GPUS=${2:-0,1,2,3}; HERE=$(cd "$(dirname "$0")" && pwd)
command -v boltz >/dev/null || { echo HOLD boltz missing; exit 2; }
CACHE="$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache"; [[ -d "$CACHE" ]] && export BOLTZ_CACHE="$CACHE"
python "$HERE/scripts/run_boltz_branch.py" --manifest "$HERE/outputs/ARCH5_STRUCTURE_MANIFEST.tsv" --structures-dir "$HERE/outputs/arch5_structures" --work-dir "$HERE/work/arch5_boltz" --logs-dir "$HERE/logs/arch5_boltz" --gpus "$GPUS" --status "$HERE/outputs/ARCH5_BOLTZ_STATUS.tsv"
echo PASS_STAGE04B9B_ARCH5_BOLTZ5_AVAILABLE
