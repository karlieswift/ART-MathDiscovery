#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
python "$HERE/tests/run_selftest.py"
python -m py_compile "$HERE"/scripts/*.py
bash -n "$HERE"/*.sh
