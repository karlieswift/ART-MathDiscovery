#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root required}
THREADS=${2:-32}
H=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$H/logs" "$H/work" "$H/outputs"
"$H/run_selftest.sh"
python "$H/scripts/stage04b6.py" prepare --project-root "$P" --home "$H" --outdir "$H/outputs" --threads "$THREADS" 2>&1 | tee "$H/logs/prepare.log"
