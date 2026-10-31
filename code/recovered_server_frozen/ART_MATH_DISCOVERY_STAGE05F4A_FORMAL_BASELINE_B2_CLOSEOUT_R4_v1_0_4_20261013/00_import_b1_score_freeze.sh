#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ARC=${1:-}
if [[ -z "$ARC" ]]; then echo "Usage: $0 /absolute/path/to/ART_MATH_DISCOVERY_STAGE05F4A_B1_R2_SCORE_FREEZE_*.tar.gz" >&2; exit 2; fi
python "$HERE/scripts/import_b1.py" --archive "$ARC" | tee "$HERE/logs/run_00_import_b1.log"
