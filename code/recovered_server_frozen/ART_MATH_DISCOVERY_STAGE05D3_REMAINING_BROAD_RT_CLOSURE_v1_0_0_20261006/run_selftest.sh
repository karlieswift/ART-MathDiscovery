#!/usr/bin/env bash
set -euo pipefail
python scripts/selftest.py "${1:-}"
