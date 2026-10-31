#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo "usage: $0 <E3D-results-tar-or-dir> <project-root> [cpu]" >&2
  exit 2
fi
E3D="$1"; P="$2"; CPU="${3:-16}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
rm -rf outputs work
mkdir -p outputs work logs
if ! command -v phmmer >/dev/null 2>&1; then
  echo "[HOLD] phmmer not found. Activate the project art_rt_evidence environment (HMMER 3.4)." >&2
  exit 3
fi
python scripts/stage05e3e.py --e3d "$E3D" --project-root "$P" --cpu "$CPU" 2>&1 | tee logs/stage05e3e.log
