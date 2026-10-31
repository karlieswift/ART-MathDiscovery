#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 3 ]; then echo "usage: $0 <E3C_results.tar.gz|dir> <E3A_results.tar.gz|dir> <project_root>" >&2; exit 2; fi
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
rm -rf outputs work; mkdir -p outputs work logs
set +e
python scripts/stage05e3d.py --e3c "$1" --e3a "$2" --project-root "$3" --outdir outputs --workdir work 2>&1 | tee logs/stage05e3d.log
rc=${PIPESTATUS[0]}
set -e
if [ "$rc" -ne 0 ]; then echo "[HOLD] Stage05E3D v1.0.5 exited rc=$rc" >&2; exit "$rc"; fi
cat outputs/STAGE05E3D_GATE.json
