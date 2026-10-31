#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
JOBS="${1:-4}"
cd "$ROOT"
mkdir -p outputs logs work
exec 9>work/.run.lock
flock -n 9 || { echo '[HOLD] Another run holds the write lock.'; exit 2; }
python scripts/verify.py --root "$ROOT"
# Source retrieval may leave an explicit HOLD. It does not prevent useful RNA
# characterization and it never authorizes a catalogue or novelty PASS.
python scripts/author_sources.py --root "$ROOT" 2>&1 | tee logs/author_sources.log
set +e
python scripts/fold.py --root "$ROOT" --jobs "$JOBS" 2>&1 | tee logs/RNAfold.log
RC=${PIPESTATUS[0]}
set -e
python scripts/finish.py --root "$ROOT"
exit "$RC"
