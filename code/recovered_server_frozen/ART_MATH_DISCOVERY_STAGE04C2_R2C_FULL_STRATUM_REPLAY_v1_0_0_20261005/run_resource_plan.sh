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
"$PY" scripts/run_replay.py plan --stage04a-root "$A04" --ug27-root "$UG"
