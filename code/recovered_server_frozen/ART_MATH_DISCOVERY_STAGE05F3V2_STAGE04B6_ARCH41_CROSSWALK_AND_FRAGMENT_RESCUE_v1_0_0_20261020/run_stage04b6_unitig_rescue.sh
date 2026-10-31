#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root required}
WORKERS=${2:-4}
H=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$H/logs" "$H/unitigs" "$H/outputs" "$H/work"
test -f "$H/outputs/STAGE04B6_PREPARE_GATE.json" || { echo HOLD_RUN_PREPARE_FIRST; exit 2; }
python "$H/scripts/stage04b6.py" download --home "$H" --outdir "$H/outputs" --unitig-dir "$H/unitigs" --workers "$WORKERS" 2>&1 | tee "$H/logs/unitig_download.log"
python "$H/scripts/stage04b6.py" rescue --project-root "$P" --home "$H" --outdir "$H/outputs" --unitig-dir "$H/unitigs" 2>&1 | tee "$H/logs/unitig_rescue.log"
python "$H/scripts/stage04b6.py" finalize --home "$H" --outdir "$H/outputs" 2>&1 | tee "$H/logs/finalize.log"
