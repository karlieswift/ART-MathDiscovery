#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
DOWNLOAD=0
if [ "${1:-}" = "--download" ]; then DOWNLOAD=1; fi
echo '[0/3] registry validation'
python scripts/validate_registry.py | tee logs/00_registry_validation.log
if [ "$DOWNLOAD" -eq 1 ]; then
  echo '[1/3] acquire currently resolved/candidate sequence records'
  scripts/acquire_sequences.sh 2>&1 | tee logs/01_acquire_sequences.log
else
  echo '[1/3] acquisition skipped (registry-only freeze)'
fi
echo '[2/3] final gate'
python scripts/finalize.py | tee logs/02_finalize.log
echo '[PASS] Stage04C0 registry freeze complete'
