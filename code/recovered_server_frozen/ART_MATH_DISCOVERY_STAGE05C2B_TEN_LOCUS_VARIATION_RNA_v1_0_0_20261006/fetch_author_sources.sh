#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p outputs logs work
exec 9>work/.run.lock
flock -n 9 || { echo '[HOLD] Another run holds the write lock.'; exit 2; }
python scripts/verify.py --root "$ROOT"
ARGS=()
if [ "$#" -gt 0 ]; then ARGS=(--pdf "$1"); fi
python scripts/author_sources.py --root "$ROOT" "${ARGS[@]}" 2>&1 | tee logs/author_sources.log
