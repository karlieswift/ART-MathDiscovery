#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_MATHDISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
R=${1:-${F3V2_SOURCE_R1_STAGE:-$P/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_R1_v1_0_1_20261012}}
python "$S/scripts/merge_and_freeze_r2.py" --stage "$S" --source-r1-stage "$R"
