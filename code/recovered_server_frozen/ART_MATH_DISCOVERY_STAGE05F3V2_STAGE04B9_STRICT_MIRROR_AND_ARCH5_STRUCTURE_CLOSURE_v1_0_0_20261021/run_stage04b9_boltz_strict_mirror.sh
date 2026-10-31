#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root}; GPUS=${2:-0,1,2,3}; HERE=$(cd "$(dirname "$0")" && pwd)
command -v boltz >/dev/null || { echo HOLD boltz missing; exit 2; }
CACHE="$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache"; [[ -d "$CACHE" ]] && export BOLTZ_CACHE="$CACHE"
python "$HERE/scripts/run_boltz_branch.py" --manifest "$HERE/outputs/STRICT_MIRROR_STRUCTURE_MANIFEST_32.tsv" --structures-dir "$HERE/outputs/strict_mirror_structures" --work-dir "$HERE/work/strict_mirror_boltz" --logs-dir "$HERE/logs/strict_mirror_boltz" --gpus "$GPUS" --status "$HERE/outputs/STRICT_MIRROR_BOLTZ_STATUS.tsv"
echo PASS_STAGE04B9A_STRICT_MIRROR_BOLTZ32_AVAILABLE
