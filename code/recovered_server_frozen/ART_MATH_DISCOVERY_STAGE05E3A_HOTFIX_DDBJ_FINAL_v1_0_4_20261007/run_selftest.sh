#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python - "$ROOT" <<'PY'
import sys
from pathlib import Path
p=Path(sys.argv[1])/'scripts/recover_ok905446_ddbj_v104.py'
s=p.read_text()
assert 'TARGET = "OK905446.1"' in s
assert 'database=na' in s
assert 'format=flatfile' in s
assert 'filetype=text' in s
assert 'DDBJ_GETENTRY_INSDC' in s
assert 'scientific_settings_changed' in s
print('[PASS] Stage05E3A v1.0.4 DDBJ final-record hotfix selftest')
PY
