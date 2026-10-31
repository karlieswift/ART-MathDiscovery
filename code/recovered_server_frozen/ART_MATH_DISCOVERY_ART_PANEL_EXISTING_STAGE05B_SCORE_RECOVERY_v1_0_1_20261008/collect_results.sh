#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
test -f outputs/ART_PANEL_STAGEA1_GATE.json || { echo "missing gate" >&2; exit 2; }
python - <<'PY'
import json
x=json.load(open('outputs/ART_PANEL_STAGEA1_GATE.json'))
assert x.get('gate_pass') is True, x
PY
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_PANEL_STAGEA1_EXISTING_STAGE05B_SCORE_RECOVERY_RESULTS_${TS}.tar.gz"
tar -czf "$PKG" outputs logs config docs scripts README.md run_recover_existing_scores.sh collect_results.sh tests
sha256sum "$PKG" > "$PKG.sha256"
echo "$HERE/$PKG"
cat "$PKG.sha256"
