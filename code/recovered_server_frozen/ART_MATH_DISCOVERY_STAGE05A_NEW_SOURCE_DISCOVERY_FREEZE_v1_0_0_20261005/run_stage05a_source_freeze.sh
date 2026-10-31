#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
command -v datasets >/dev/null || { echo '[STOP] NCBI datasets is required; use the drhu environment on r750-1 where it is already available.'; exit 2; }
datasets version
mkdir -p logs
"$PY" scripts/discovery_source.py acquire "$@" 2>&1 | tee logs/source_acquisition.log
