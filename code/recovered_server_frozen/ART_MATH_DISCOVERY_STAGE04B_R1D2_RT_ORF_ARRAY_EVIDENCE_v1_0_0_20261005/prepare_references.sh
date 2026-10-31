#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
exec 9> .execution.lock
flock -n 9 || { echo '[HOLD] another process owns this run directory'; exit 2; }
"${PYTHON:-python3}" scripts/verify_package.py
"${PYTHON:-python3}" -u scripts/prepare_refs.py "$@"
