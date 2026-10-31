#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
"$ROOT/scripts/run_sensitive_homology.sh" "$@" 2>&1 | tee "$ROOT/logs/run_stage05d4r2.log"
