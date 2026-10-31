#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
python - "$HERE/scripts/stage04b8r1.py" <<'PY'
import ast,sys
p=sys.argv[1]
ast.parse(open(p).read())
print("PASS_STAGE04B8R1_PYTHON_COMPILE")
PY
bash -n "$HERE/run_stage04b8r1_audit.sh"
echo PASS_STAGE04B8R1_SELFTEST
