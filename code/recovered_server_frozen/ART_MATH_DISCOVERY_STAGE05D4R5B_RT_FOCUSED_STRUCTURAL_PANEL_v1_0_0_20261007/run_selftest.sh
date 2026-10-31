#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
python -m py_compile scripts/build_rt_panel.py scripts/summarize_r5b.py
python tests/selftest.py
