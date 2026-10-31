#!/usr/bin/env bash
set -euo pipefail
E3A="${1:?usage: run_stage05e3a_fetchfix_v101.sh <E3A_ROOT>}"
cd "$E3A"

# Preserve v1.0.0 invalid transport payloads before retrying.
mkdir -p raw/genbank_v100_invalid raw/genbank outputs logs
for f in raw/genbank/*.gbff; do
  [ -e "$f" ] || continue
  if ! grep -aq '^LOCUS' "$f"; then
    mv -f "$f" "raw/genbank_v100_invalid/$(basename "$f")"
  fi
done

python "$(cd "$(dirname "$0")" && pwd)/fetch_genbank_v101.py" \
  --ids frozen/STAGE05E3A_SEQUENCE_IDS.txt \
  --outdir raw/genbank \
  --manifest raw/STAGE05E3A_GENBANK_RETRIEVAL_MANIFEST_v101.tsv \
  --diagnostics raw/STAGE05E3A_GENBANK_RETRIEVAL_DIAGNOSTICS_v101.tsv \
  2>&1 | tee logs/fetch_genbank_v101.log

python scripts/analyze_official_context.py \
  --frozen frozen \
  --gbdir raw/genbank \
  --outdir outputs \
  2>&1 | tee logs/analyze_official_context_v101.log

python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/GATE_STAGE05E3A.json')
if not p.exists(): raise SystemExit('Missing Stage05E3A gate after v1.0.1 retrieval hotfix')
x=json.loads(p.read_text())
if not x.get('gate_pass'): raise SystemExit('Stage05E3A gate did not pass: '+repr(x))
print('[PASS] Stage05E3A v1.0.1 retrieval hotfix -> official-context analysis complete')
print(json.dumps(x,indent=2))
PY
