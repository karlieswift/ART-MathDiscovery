#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
E3F="${1:?E3F result tar required}"
P="${2:?project root required}"
rm -rf "$HERE/outputs" "$HERE/work"
mkdir -p "$HERE/outputs" "$HERE/work"
python "$HERE/scripts/stage05e3g.py" --e3f "$E3F" --project-root "$P" --out "$HERE/outputs" --work "$HERE/work"
python - <<'PY' "$HERE/outputs/STAGE05E3G_GATE.json"
import json,sys
j=json.load(open(sys.argv[1]))
assert j['gate_pass'],j
print('[PASS] Stage05E3G orthogonal evidence recovery gate')
PY
