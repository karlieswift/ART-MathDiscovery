#!/usr/bin/env bash
set -euo pipefail
P="${1:?project root required}"
THREADS="${2:-16}"
mkdir -p outputs logs work
python scripts/stage05f1.py "$P" --threads "$THREADS" --outdir outputs 2>&1 | tee logs/stage05f1.log
