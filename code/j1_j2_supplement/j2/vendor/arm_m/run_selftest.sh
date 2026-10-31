#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
python -m py_compile scripts/*.py tests/*.py
python tests/test_engine.py
python scripts/00_freeze_engine.py --root "$ROOT" >/tmp/stage03a_engine_freeze_selftest.json
python - <<'PY'
import json
x=json.load(open('outputs/STAGE03A_ENGINE_FREEZE_GATE.json'))
assert x['gate_pass'] and x['gate']=='PASS_STAGE03A_FROZEN_ENGINE_BEFORE_HOLDOUT_SCORING'
print('[PASS] Stage03A engine-freeze gate selftest')
PY
rm -f outputs/STAGE03A_ENGINE_FREEZE_GATE.json
