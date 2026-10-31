#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
python - <<'PY' "$HERE/outputs/STAGE05E3H_GATE.json"
import json,sys
x=json.load(open(sys.argv[1])); assert x['gate_pass'] is True
PY
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="$HERE/ART_STAGE05E3H_GLOBAL_NOVELTY_SHORT_PROTEIN_CLOSURE_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" -C "$HERE" README.md config docs scripts tests run_stage05e3h.sh collect_stage05e3h.sh outputs logs
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$OUT"
