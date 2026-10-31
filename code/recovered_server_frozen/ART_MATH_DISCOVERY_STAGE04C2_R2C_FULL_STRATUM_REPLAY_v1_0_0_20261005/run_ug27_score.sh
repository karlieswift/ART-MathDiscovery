#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
PY="${PYTHON:-python}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p outputs logs
A04="${1:?provide the Stage04A root}"
UG="${2:?provide the UG27 source-freeze root}"
THREADS="${3:-8}"
"$PY" scripts/run_replay.py score --stage04a-root "$A04" --ug27-root "$UG" --threads "$THREADS" 2>&1 | tee -a logs/ug27_run.log
