#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
STR=${1:-}
[[ "$STR" =~ ^[ABCD]$ ]] || { echo "usage: $0 A|B|C|D [SOURCE_R1_STAGE]" >&2; exit 2; }
P=${ART_MATHDISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
R=${2:-${F3V2_SOURCE_R1_STAGE:-$P/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_R1_v1_0_1_20261012}}
mkdir -p "$S/logs"
python "$S/scripts/acquire_one_stratum_r2.py" --stage "$S" --source-r1-stage "$R" --stratum "$STR"
