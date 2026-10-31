#!/usr/bin/env bash
set -euo pipefail
python -m py_compile "$(cd "$(dirname "$0")/.." && pwd)/scripts/weak_family_repeat_control.py"
echo '[PASS] F3B4B1C syntax selftest'
