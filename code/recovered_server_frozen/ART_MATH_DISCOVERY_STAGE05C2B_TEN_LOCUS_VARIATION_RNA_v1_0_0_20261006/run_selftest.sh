#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
PYTHONPATH="$ROOT/scripts${PYTHONPATH:+:$PYTHONPATH}" python -m unittest discover -s tests -v
python scripts/verify.py --root "$ROOT"
