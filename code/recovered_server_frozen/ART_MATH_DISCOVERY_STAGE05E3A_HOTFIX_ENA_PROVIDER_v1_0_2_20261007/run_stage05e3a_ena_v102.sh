#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then
  echo "usage: $0 <E3A_ROOT>" >&2
  exit 2
fi
E3A="$(cd "$1" && pwd)"
FIX="$(cd "$(dirname "$0")" && pwd)"
cd "$E3A"

./run_selftest.sh

# Preserve blocked/invalid NCBI payloads and manifests for provenance.
mkdir -p raw/genbank_ncbi_blocked_v100_v101 raw/genbank outputs logs
if find raw/genbank -maxdepth 1 -type f | grep -q .; then
  find raw/genbank -maxdepth 1 -type f -exec mv -t raw/genbank_ncbi_blocked_v100_v101 {} + || true
fi
for f in raw/STAGE05E3A_GENBANK_RETRIEVAL_MANIFEST.tsv raw/STAGE05E3A_GENBANK_RETRIEVAL_DIAGNOSTICS_v101.tsv; do
  if [ -f "$f" ]; then
    mv "$f" raw/genbank_ncbi_blocked_v100_v101/ || true
  fi
done

rm -f raw/genbank/*.gbff

python "$FIX/scripts/fetch_ena_to_genbank_v102.py" \
  --ids frozen/STAGE05E3A_SEQUENCE_IDS.txt \
  --outdir raw/genbank \
  --manifest raw/STAGE05E3A_ENA_RETRIEVAL_MANIFEST_v102.tsv \
  --diagnostics raw/STAGE05E3A_ENA_RETRIEVAL_DIAGNOSTICS_v102.tsv \
  2>&1 | tee logs/fetch_ena_v102.log

N=$(find raw/genbank -maxdepth 1 -type f -name '*.gbff' | wc -l | tr -d ' ')
if [ "$N" != "61" ]; then
  echo "[FAIL] expected 61 validated GenBank flatfiles after ENA conversion; got $N" >&2
  exit 1
fi

echo "[PASS] ENA retrieval/validation/conversion complete: 61/61"

# Downstream biological analysis is unchanged from Stage05E3A v1.0.0.
python scripts/analyze_official_context.py \
  --frozen frozen \
  --gbdir raw/genbank \
  --outdir outputs \
  2>&1 | tee logs/analyze_official_context_v102.log

python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/GATE_STAGE05E3A.json')
if not p.exists(): raise SystemExit('[FAIL] missing E3A gate')
x=json.loads(p.read_text())
if not x.get('gate_pass'): raise SystemExit(f"[FAIL] E3A gate: {x}")
if x.get('n_genbank_records_parsed') != 61: raise SystemExit(f"[FAIL] parsed {x.get('n_genbank_records_parsed')} records, expected 61")
print('[PASS] Stage05E3A official-context analysis complete after ENA provider hotfix')
PY
