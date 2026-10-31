#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
"${PYTHON:-python3}" scripts/verify_package.py
"${PYTHON:-python3}" -m unittest discover -s tests -v
