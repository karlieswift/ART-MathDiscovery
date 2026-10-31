#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then
  echo "usage: $0 /absolute/path/to/Stage04A_root" >&2
  exit 2
fi
ROOT=$(cd "$(dirname "$0")" && pwd)
python "$ROOT/scripts/run_stage04b_r1.py" "$1" 2>&1 | tee "$ROOT/logs/STAGE04B_R1_RUN.log"
