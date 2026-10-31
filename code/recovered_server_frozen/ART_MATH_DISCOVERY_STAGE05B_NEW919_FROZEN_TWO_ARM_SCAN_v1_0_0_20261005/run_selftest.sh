#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python -c 'import sys; assert sys.version_info >= (3,9), "Python >=3.9 required; use bioh200"; import numpy, pandas'
python -m unittest discover -s tests -p 'test_*.py' -v
