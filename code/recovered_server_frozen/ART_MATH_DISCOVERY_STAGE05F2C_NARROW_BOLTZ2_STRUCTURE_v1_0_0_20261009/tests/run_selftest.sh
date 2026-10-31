#!/usr/bin/env bash
set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
p=Path('config/FROZEN_TARGET.json')
d=json.loads(p.read_text())
assert d['target_query_id']=='SF_00172800'
assert d['target_aa_length']==205
assert d['target_rvt_profile']=='RVT-UG9'
assert len(d['expected_accessions'])==2
assert d['boltz_protocol']['seed']==20261007
print('[PASS] Stage05F2C frozen-target/protocol selftest')
PY
