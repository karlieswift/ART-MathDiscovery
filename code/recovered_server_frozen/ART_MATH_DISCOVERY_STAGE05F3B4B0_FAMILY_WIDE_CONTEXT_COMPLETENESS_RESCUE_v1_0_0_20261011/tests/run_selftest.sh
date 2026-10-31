#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
python - "$HERE/scripts/family_wide_context_rescue.py" <<'PY'
import ast,sys
p=sys.argv[1]
ast.parse(open(p).read())
print('[PASS] F3B4B0 python syntax selftest')
PY
