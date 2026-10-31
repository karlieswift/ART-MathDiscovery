#!/usr/bin/env bash
set -euo pipefail
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
S=$(cd "$(dirname "$0")" && pwd)
B4D="$P/ART_MATH_DISCOVERY_STAGE05F3B4D_SEVEN_FAMILY_A_TO_E_CLOSEOUT_v1_0_0_20261011"
mkdir -p "$S"/{outputs,inputs,frozen,logs,work,tmp}
python "$S/scripts/prepare.py" "$B4D" "$S/outputs" "$S/inputs" "$S/frozen" "$S/frozen/F3B4E_IMPLEMENTATION_LOCK.json" | tee "$S/logs/prepare.log"
echo PASS_STAGE05F3B4E_PREPARE_RUN_COMPLETE
