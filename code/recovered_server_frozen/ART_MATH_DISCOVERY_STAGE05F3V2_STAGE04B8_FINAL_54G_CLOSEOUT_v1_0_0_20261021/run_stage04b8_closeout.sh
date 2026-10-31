#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
P="${1:?project root required}"
mkdir -p "$HERE/logs" "$HERE/outputs"
"$HERE/run_selftest.sh"
python "$HERE/scripts/stage04b8.py" closeout --project-root "$P" --home "$HERE" \
  2>&1 | tee "$HERE/logs/stage04b8_closeout.log"
