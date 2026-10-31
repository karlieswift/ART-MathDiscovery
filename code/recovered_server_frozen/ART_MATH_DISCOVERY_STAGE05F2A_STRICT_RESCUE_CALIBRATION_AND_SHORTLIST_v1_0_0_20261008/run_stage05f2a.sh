#!/usr/bin/env bash
set -euo pipefail
F1=${1:?Stage05F1 directory or results tar required}
R1=${2:?Stage05F1R1 directory or results tar required}
cd "$(dirname "$0")"
rm -rf outputs work logs
mkdir -p outputs work logs
python scripts/stage05f2a.py "$F1" "$R1" 2>&1 | tee logs/stage05f2a.log
