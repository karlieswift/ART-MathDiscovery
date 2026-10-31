#!/usr/bin/env bash
set -euo pipefail
LOG=${1:-logs/boltz_p1.log}
[ -f "$LOG" ] || { echo "missing $LOG"; exit 2; }
echo '=== first fatal/error context ==='
grep -n -m 20 -E 'Traceback|RuntimeError|CUDA out of memory|OutOfMemory|ChildFailedError|ProcessRaisedException|SIGTERM|SIGKILL|NCCL|Broken pipe|failed examples' "$LOG" || true
echo '=== tail ==='
tail -120 "$LOG"
