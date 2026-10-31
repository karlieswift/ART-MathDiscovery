#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
PY="${PYTHON:-python}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p outputs logs
sha256sum -c PACKAGE_MANIFEST.sha256
"$PY" -m unittest discover -s tests -p "test_*.py" -v 2>&1 | tee logs/selftest.log
