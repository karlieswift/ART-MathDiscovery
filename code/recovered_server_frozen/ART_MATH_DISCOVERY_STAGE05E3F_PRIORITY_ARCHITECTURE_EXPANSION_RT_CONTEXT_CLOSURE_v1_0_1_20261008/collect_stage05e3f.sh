#!/usr/bin/env bash
set -euo pipefail
python - <<'PY'
import json
p='outputs/STAGE05E3F_GATE.json'
g=json.load(open(p))
assert g['gate_pass'] and g['gate']=='PASS_E3F_PRIORITY_ARCHITECTURE_EXPANSION_AND_RT_CONTEXT_CLOSURE'
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05E3F_PRIORITY_ARCHITECTURE_EXPANSION_RT_CONTEXT_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" README.md config docs scripts tests run_stage05e3f.sh collect_stage05e3f.sh outputs
sha256sum "$OUT" > "$OUT.sha256"
echo "$PWD/$OUT"
cat "$OUT.sha256"
