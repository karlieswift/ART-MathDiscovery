#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 3 ]; then
  echo "usage: $0 <E3C_results.tar.gz|dir> <E3A_results.tar.gz|dir> <project_root>" >&2
  exit 2
fi
E3C="$1"; E3A="$2"; P="$3"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
rm -rf outputs work
mkdir -p outputs work logs
set +e
python scripts/stage05e3d.py --e3c "$E3C" --e3a "$E3A" --project-root "$P" --outdir outputs --workdir work > logs/stage05e3d.log 2>&1
rc=$?
set -e
cat logs/stage05e3d.log
if [ "$rc" -ne 0 ]; then
  echo "[HOLD] Stage05E3D exited rc=$rc; inspect outputs/STAGE05E3D_SCHEMA_DIAGNOSTICS.tsv and summary." >&2
  exit "$rc"
fi
echo "[PASS] Stage05E3D exact sequence materialization gate"
