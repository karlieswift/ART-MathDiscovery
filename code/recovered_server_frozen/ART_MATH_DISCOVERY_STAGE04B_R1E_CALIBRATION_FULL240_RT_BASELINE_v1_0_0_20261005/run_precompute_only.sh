#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python scripts/run_stage.py --precompute-only
