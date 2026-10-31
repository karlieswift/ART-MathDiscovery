#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then
  echo "usage: $0 <E3A_ROOT>" >&2
  exit 2
fi
E3A="$(cd "$1" && pwd)"
FIX="$(cd "$(dirname "$0")" && pwd)"
cd "$E3A"

python - <<'PY'
import csv, json
from pathlib import Path

gp=Path('outputs/GATE_STAGE05E3A.json')
sp=Path('outputs/STAGE05E3A_SUMMARY.json')
mp=Path('raw/STAGE05E3A_COMBINED_PROVIDER_MANIFEST_v104.tsv')
for p in (gp,sp,mp):
    if not p.exists():
        raise SystemExit(f'[FAIL] required Stage05E3A closeout file missing: {p}')
g=json.loads(gp.read_text())
s=json.loads(sp.read_text())
assert g.get('gate_pass') is True, g
assert g.get('gate') == 'PASS_E2A_OFFICIAL_CONTEXT_AND_LINEAGE_AUDIT_COMPLETE', g
assert int(g.get('n_genbank_records_parsed', -1)) == 61, g
assert int(s.get('n_genbank_records_parsed', -1)) == 61, s
with mp.open() as h:
    rows=list(csv.DictReader(h, delimiter='\t'))
if len(rows) != 61:
    raise SystemExit(f'[FAIL] combined provider manifest rows={len(rows)} expected=61')
if any(r.get('status') != 'VALID' for r in rows):
    bad=[r.get('sequence_id') for r in rows if r.get('status')!='VALID']
    raise SystemExit(f'[FAIL] combined provider manifest has non-VALID rows: {bad[:10]}')
ids=[r.get('sequence_id','') for r in rows]
if len(set(ids)) != 61:
    raise SystemExit('[FAIL] combined provider manifest sequence_id uniqueness != 61')
print('[PASS] Stage05E3A closeout gate + 61/61 combined provider provenance')
PY

STAMP=$(date +%Y%m%d_%H%M%S)
OUT="../ART_STAGE05E3A_E2A_OFFICIAL_CONTEXT_LINEAGE_RESULTS_R1_${STAMP}.tar.gz"
TMP_LIST="$(mktemp)"
trap 'rm -f "$TMP_LIST"' EXIT

# Core frozen scientific package and completed outputs.
for p in frozen config docs outputs logs README.md run_selftest.sh run_stage05e3a.sh collect_stage05e3a.sh scripts; do
  [ -e "$p" ] && printf '%s\n' "$p" >> "$TMP_LIST"
done

# Canonical multi-provider retrieval provenance after v1.0.4.
printf '%s\n' raw/STAGE05E3A_COMBINED_PROVIDER_MANIFEST_v104.tsv >> "$TMP_LIST"

# Provider-specific retrieval ledgers/diagnostics are audit evidence and are included when present.
for p in \
  raw/STAGE05E3A_ENA_RETRIEVAL_MANIFEST_v102.tsv \
  raw/STAGE05E3A_ENA_RETRIEVAL_DIAGNOSTICS_v102.tsv \
  raw/retrieval_v103 \
  raw/retrieval_v104; do
  [ -e "$p" ] && printf '%s\n' "$p" >> "$TMP_LIST"
done

# Include this collector hotfix itself as provenance.
mkdir -p hotfix_provenance_v105
cp -f "$FIX/README.md" hotfix_provenance_v105/README.md
cp -f "$FIX/collect_stage05e3a_v105.sh" hotfix_provenance_v105/collect_stage05e3a_v105.sh
printf '%s\n' hotfix_provenance_v105 >> "$TMP_LIST"

# Sort/deduplicate top-level inclusion paths for deterministic tar input ordering.
sort -u "$TMP_LIST" -o "$TMP_LIST"
tar -czf "$OUT" -T "$TMP_LIST"
sha256sum "$OUT" > "${OUT}.sha256"

python - "$OUT" <<'PY'
import tarfile, sys
p=sys.argv[1]
with tarfile.open(p,'r:gz') as t:
    names=set(t.getnames())
required={
 'outputs/GATE_STAGE05E3A.json',
 'outputs/STAGE05E3A_SUMMARY.json',
 'raw/STAGE05E3A_COMBINED_PROVIDER_MANIFEST_v104.tsv',
 'hotfix_provenance_v105/README.md',
 'hotfix_provenance_v105/collect_stage05e3a_v105.sh',
}
missing=sorted(required-names)
if missing:
    raise SystemExit(f'[FAIL] result archive missing required closeout members: {missing}')
print('[PASS] Stage05E3A R1 result archive contains gate, summary, canonical provider manifest, and collector-hotfix provenance')
PY

echo "$OUT"
cat "${OUT}.sha256"
