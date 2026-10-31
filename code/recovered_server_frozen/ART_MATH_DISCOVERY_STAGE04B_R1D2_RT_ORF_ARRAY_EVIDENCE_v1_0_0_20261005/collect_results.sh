#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
# Do not collect while a run writes output files.
exec 9> .execution.lock
flock -n 9 || { echo '[HOLD] analysis is still writing; collect after completion'; exit 2; }
"${PYTHON:-python3}" scripts/collect_results.py
