#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
python - "$ROOT" <<'PY'
import json,sys,tempfile
from pathlib import Path
root=Path(sys.argv[1]); p=json.loads((root/'config/HOTFIX_POLICY.json').read_text())
assert p['provider']=='DDBJ_GETENTRY_INSDC'
assert p['mutate_original_c0'] is False
assert p['run_acnv_or_recurrence'] is False
print('[PASS] Stage04C0-R1 hotfix policy selftest')
PY
