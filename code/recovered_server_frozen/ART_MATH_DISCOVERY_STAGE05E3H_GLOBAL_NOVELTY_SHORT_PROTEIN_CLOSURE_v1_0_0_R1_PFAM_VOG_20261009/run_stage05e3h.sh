#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 4 ]; then
  echo "usage: $0 E3G_RESULTS.tar.gz E3F_RESULTS.tar.gz E3D_RESULTS.tar.gz PROJECT_ROOT [THREADS]" >&2
  exit 2
fi
E3G="$1"; E3F="$2"; E3D="$3"; P="$4"; T="${5:-16}"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HERE/outputs" "$HERE/logs" "$HERE/work"
python "$HERE/scripts/stage05e3h.py" \
  --e3g "$E3G" --e3f "$E3F" --e3d "$E3D" --project-root "$P" --threads "$T" --outdir "$HERE/outputs" \
  2>&1 | tee "$HERE/logs/stage05e3h.log"
python - <<'PY' "$HERE/outputs/STAGE05E3H_GATE.json"
import json,sys
x=json.load(open(sys.argv[1])); assert x['gate_pass'] is True; print('[PASS] Stage05E3H accelerated local closure gate')
PY
