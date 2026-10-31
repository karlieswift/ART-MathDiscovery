#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
command -v datasets >/dev/null || { echo '[FAIL] datasets CLI not found'; exit 2; }
python "$ROOT/scripts/fetch_datasets_cache.py"
