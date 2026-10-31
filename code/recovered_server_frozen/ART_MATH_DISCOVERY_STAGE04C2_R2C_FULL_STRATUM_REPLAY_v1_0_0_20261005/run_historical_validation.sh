#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
PY="${PYTHON:-python}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p outputs logs
A04="${1:?provide the Stage04A root}"
THREADS="${2:-8}"
"$PY" scripts/run_replay.py validate --stage04a-root "$A04" --threads "$THREADS" 2>&1 | tee logs/historical_validation.log
