#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
python - <<'PY'
import json
from pathlib import Path
g=json.loads(Path('outputs/GATE_STAGE05E3A.json').read_text())
assert g.get('gate_pass') is True, g
print('[PASS] Stage05E3A gate')
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="../ART_STAGE05E3A_E2A_OFFICIAL_CONTEXT_LINEAGE_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" frozen config docs outputs raw/STAGE05E3A_GENBANK_RETRIEVAL_MANIFEST.tsv logs README.md run_selftest.sh run_stage05e3a.sh collect_stage05e3a.sh scripts
sha256sum "$OUT" > "${OUT}.sha256"
echo "$OUT"
cat "${OUT}.sha256"
