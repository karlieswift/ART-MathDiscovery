#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
python - <<'PY' "$HERE/outputs/STAGE05E3G_GATE.json"
import json,sys
j=json.load(open(sys.argv[1]))
if not j.get('gate_pass'): raise SystemExit('gate not PASS')
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="$HERE/ART_STAGE05E3G_PRIORITY_FAMILY_ORTHOGONAL_EVIDENCE_RECOVERY_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" -C "$HERE" README.md config docs outputs scripts tests run_stage05e3g.sh collect_stage05e3g.sh
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$OUT"
