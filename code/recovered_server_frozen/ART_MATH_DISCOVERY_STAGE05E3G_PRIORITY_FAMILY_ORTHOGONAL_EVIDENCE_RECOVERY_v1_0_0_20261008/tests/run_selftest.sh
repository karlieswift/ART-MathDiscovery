#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/e3f/outputs" "$TMP/P/ART_MATH_DISCOVERY_STAGE05D4_R2_SWISSPROT/outputs"
python - <<'PY' "$TMP/e3f/outputs" "$TMP/P/ART_MATH_DISCOVERY_STAGE05D4_R2_SWISSPROT/outputs"
from pathlib import Path
import sys,hashlib
out=Path(sys.argv[1]); d=Path(sys.argv[2])
seqs={'ARCH250_0010':'M'*56,'ARCH250_0013':'A'*109,'ARCH250_0033':'G'*398}
for i,(fam,s) in enumerate(seqs.items()):
    (out/f'{fam}.all_local_homologs.faa').write_text(f'>{fam}_W sequence_id=X{i} evalue=1e-20 witness=True\n{s}\n')
with (out/'STAGE05E3F_FAMILY_EXPANSION_SUMMARY.tsv').open('w') as f:
    f.write('family_id\te3f_architecture_expansion_class\n')
    for fam in seqs:f.write(f'{fam}\tTEST_CLASS\n')
for fam,s in seqs.items():
    h=hashlib.sha256(s.encode()).hexdigest()
    with (d/'hits.tsv').open('a') as f:f.write(f'{fam}\tAA_{h}\tSwissProt hit\n')
PY
tar -czf "$TMP/e3f.tar.gz" -C "$TMP/e3f" .
python "$HERE/scripts/stage05e3g.py" --e3f "$TMP/e3f.tar.gz" --project-root "$TMP/P" --out "$TMP/o" --work "$TMP/w" >/dev/null
python - <<'PY' "$TMP/o/STAGE05E3G_GATE.json" "$TMP/o/STAGE05E3G_EXISTING_ORTHOGONAL_EVIDENCE.tsv"
import json,sys,csv
j=json.load(open(sys.argv[1])); assert j['gate_pass']; assert j['priority_family_n']==3
r=list(csv.DictReader(open(sys.argv[2]),delimiter='\t')); assert len(r)>=3
print('[PASS] Stage05E3G evidence-recovery selftest')
PY
