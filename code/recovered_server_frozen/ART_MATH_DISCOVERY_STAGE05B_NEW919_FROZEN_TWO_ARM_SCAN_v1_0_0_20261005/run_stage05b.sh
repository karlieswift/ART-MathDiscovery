#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
SOURCE="${1:?usage: ./run_stage05b.sh STAGE05A_ROOT [GENOME_WORKERS=2] [THREADS=8]}"
GW="${2:-2}"
IT="${3:-8}"
cd "$ROOT"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p logs outputs
python scripts/stage05b.py run --source "$SOURCE" --genome-workers "$GW" --threads "$IT" 2>&1 | tee -a logs/run_stage05b.log
