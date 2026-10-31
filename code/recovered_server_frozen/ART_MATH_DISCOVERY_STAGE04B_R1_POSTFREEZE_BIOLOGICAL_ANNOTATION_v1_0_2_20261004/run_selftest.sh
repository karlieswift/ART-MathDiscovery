#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
python "$ROOT/tests/test_stage04b_r1.py"
