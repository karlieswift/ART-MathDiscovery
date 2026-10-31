#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
"$ROOT/run_selftest.sh"
"$ROOT/scripts/reuse_r2_records.sh" "$1"
