#!/usr/bin/env bash
set -euo pipefail
P="${1:?project root required}"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HERE/outputs" "$HERE/logs"
python "$HERE/scripts/stage04b8r1.py" "$P" --outdir "$HERE/outputs" \
  2>&1 | tee "$HERE/logs/stage04b8r1.log"
