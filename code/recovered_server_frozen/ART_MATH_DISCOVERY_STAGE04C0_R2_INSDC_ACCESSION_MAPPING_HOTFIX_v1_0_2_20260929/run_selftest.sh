#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
python "$ROOT/tests/test_policy.py"
python "$ROOT/tests/test_validator_synthetic.py"
