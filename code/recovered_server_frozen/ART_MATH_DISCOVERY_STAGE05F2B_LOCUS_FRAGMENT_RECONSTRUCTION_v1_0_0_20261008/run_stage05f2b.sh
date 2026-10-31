#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 3 ]; then
  echo "usage: $0 <Stage05F1 live dir or result tar> <Stage05F2A live dir or result tar> <project_root> [threads]" >&2
  exit 2
fi
F1="$1"; F2A="$2"; P="$3"; T="${4:-16}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
rm -rf outputs work logs
mkdir -p outputs work logs
python scripts/stage05f2b.py "$F1" "$F2A" "$P" --threads "$T" --outdir outputs 2>&1 | tee logs/stage05f2b.log
