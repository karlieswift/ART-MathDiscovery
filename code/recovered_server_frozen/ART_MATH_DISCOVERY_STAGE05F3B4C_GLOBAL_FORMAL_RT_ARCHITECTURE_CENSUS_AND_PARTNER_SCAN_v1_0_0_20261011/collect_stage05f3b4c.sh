#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
[[ -s "$HERE/outputs/STAGE05F3B4C1_GATE.json" ]] || { echo '[FATAL] C1 gate missing'; exit 2; }
[[ -s "$HERE/outputs/STAGE05F3B4C2_GATE.json" ]] || { echo '[FATAL] C2 gate missing'; exit 2; }
python - "$HERE/outputs/STAGE05F3B4C1_GATE.json" "$HERE/outputs/STAGE05F3B4C2_GATE.json" <<'PY'
import json,sys
for p in sys.argv[1:]:
 g=json.load(open(p)); assert g['gate_pass'] is True,(p,g)
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="$HERE/../ART_STAGE05F3B4C_GLOBAL_FORMAL_RT_ARCHITECTURE_CENSUS_AND_PARTNER_SCAN_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" -C "$HERE" README_CN.txt frozen outputs logs
sha256sum "$OUT" > "$OUT.sha256"
echo "$OUT"
cat "$OUT.sha256"
