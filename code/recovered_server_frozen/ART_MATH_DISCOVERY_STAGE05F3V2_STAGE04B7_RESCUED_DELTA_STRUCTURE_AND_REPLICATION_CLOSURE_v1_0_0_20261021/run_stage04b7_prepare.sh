#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root}; THREADS=${2:-32}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/logs" "$HERE/outputs" "$HERE/work" "$HERE/inputs"
"$HERE/run_selftest.sh"
python "$HERE/scripts/stage04b7.py" prepare --project-root "$P" --home "$HERE" --threads "$THREADS"
