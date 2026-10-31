#!/usr/bin/env bash
set -euo pipefail
CPU=${1:-16}
HERE=$(cd "$(dirname "$0")" && pwd)
python "$HERE/scripts/run_score_freeze.py" --cpu "$CPU" | tee "$HERE/logs/run_01_score_freeze.log"
