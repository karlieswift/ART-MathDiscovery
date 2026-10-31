#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
PY="${PYTHON:-python}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p outputs logs
STAMP="$(date +%Y%m%d_%H%M%S)"
PKG="${ROOT##*/}_RESULTS_${STAMP}.tar.gz"
"$PY" scripts/collect_results.py "$ROOT" "$PKG"
sha256sum "$PKG" > "$PKG.sha256"
sha256sum -c "$PKG.sha256"
ls -lh "$PKG" "$PKG.sha256"
