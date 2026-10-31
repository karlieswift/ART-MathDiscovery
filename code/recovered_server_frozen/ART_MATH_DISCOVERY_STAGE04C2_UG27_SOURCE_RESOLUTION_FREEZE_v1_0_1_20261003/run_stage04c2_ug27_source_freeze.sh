#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$ROOT/outputs" "$ROOT/logs"
"$ROOT/run_selftest.sh" | tee "$ROOT/logs/selftest.log"
set +e
python "$ROOT/scripts/resolve_ug27_sources.py" --root "$ROOT" 2>&1 | tee "$ROOT/logs/resolve_ug27_sources.log"
RC=${PIPESTATUS[0]}
set -e
# RC=4 is an intentional scientific HOLD, not a software crash.
if [[ "$RC" -ne 0 && "$RC" -ne 4 ]]; then exit "$RC"; fi
cat "$ROOT/outputs/GATE.json"
exit 0
