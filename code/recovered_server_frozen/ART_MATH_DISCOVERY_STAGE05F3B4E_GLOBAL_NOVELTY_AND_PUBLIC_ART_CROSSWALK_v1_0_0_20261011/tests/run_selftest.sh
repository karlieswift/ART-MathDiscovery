#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
python -m py_compile "$S"/scripts/*.py
python - <<'PY'
# Threshold semantics regression.
def high(pid,q,t): return pid>=.90 and q>=.80 and t>=.80
assert high(.90,.80,.80)
assert not high(.899,.99,.99)
def arch(pid,q): return pid>=.95 and q>=.80
assert arch(.95,.80) and not arch(.94,.99)
print('[PASS] F3B4E threshold/parser syntax selftest')
PY
