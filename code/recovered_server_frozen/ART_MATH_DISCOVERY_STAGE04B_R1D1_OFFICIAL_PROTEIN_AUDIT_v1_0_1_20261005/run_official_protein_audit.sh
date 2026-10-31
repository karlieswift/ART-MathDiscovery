#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
PY="${PYTHON:-python}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p outputs logs
PROJECT="${1:?provide ART_MathDiscovery_v1 project root}"
shift
"$PY" scripts/audit_proteins.py --project-root "$PROJECT" "$@" 2>&1 | tee -a logs/protein_audit.log
