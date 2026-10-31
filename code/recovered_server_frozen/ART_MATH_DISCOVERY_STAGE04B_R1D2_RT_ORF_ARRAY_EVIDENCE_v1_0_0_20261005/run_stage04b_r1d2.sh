#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
CPU="${1:-4}"
cd "$ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p logs
# flock prevents two writers using the same work/output directories.
exec 9> .execution.lock
flock -n 9 || { echo '[HOLD] another process owns this run directory'; exit 2; }
"${PYTHON:-python3}" scripts/verify_package.py
"${PYTHON:-python3}" -u scripts/run_evidence.py --cpus "$CPU" 2>&1 | tee -a logs/run_stage04b_r1d2.log
