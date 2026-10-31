#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
PY="${PYTHON:-python}"
"$PY" -c 'import sys; assert sys.version_info >= (3,9), "Python >=3.9 required"'
"$PY" -m unittest discover -s tests -p 'test_*.py' -v
