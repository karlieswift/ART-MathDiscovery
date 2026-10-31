#!/usr/bin/env bash
set -euo pipefail
P=${1:?usage: bash run_stage04b1_delta_replication.sh PROJECT_ROOT}
HERE=$(cd "$(dirname "$0")" && pwd)
A="$P/ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04A_60G_SATURATION_CENSUS_v1_0_0_20261016"
D="$P/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_DISTRIBUTED_ACQUISITION_R2_v1_0_2_20261013"
OUT="$HERE/outputs"
mkdir -p "$OUT"
python "$HERE/scripts/delta_family_replication_census.py" \
  --stage04a-root "$A" \
  --delta-root "$D" \
  --out "$OUT"
