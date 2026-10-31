#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then echo "usage: $0 <E3A_ROOT>" >&2; exit 2; fi
E3A="$(cd "$1" && pwd)"
FIX="$(cd "$(dirname "$0")" && pwd)"
cd "$E3A"
./run_selftest.sh
mkdir -p raw/genbank raw/retrieval_v103 outputs logs

python - "$E3A" <<'PY'
import sys, json
from pathlib import Path
from Bio import SeqIO
E=Path(sys.argv[1]); ids=[x.strip() for x in (E/'frozen/STAGE05E3A_SEQUENCE_IDS.txt').read_text().splitlines() if x.strip()]
exp=set('''NC_005083.2
NC_019405.1
NC_019406.1
NC_019407.1
NC_019408.1
NC_019410.1
NC_023568.1
NC_028829.1
NC_031039.1
NC_031049.1
NC_043027.1
NC_048045.1
NC_048046.1
NC_048047.1
NC_048048.2
NC_052983.1
OK905446.1'''.split())
valid=[]; missing=[]
for sid in ids:
 p=E/'raw/genbank'/f'{sid}.gbff'
 ok=False
 if p.exists():
  try:
   rs=list(SeqIO.parse(str(p),'genbank')); ok=(len(rs)==1 and len(rs[0].seq)>0)
  except Exception: ok=False
 if ok: valid.append(sid)
 else: missing.append(sid)
print(json.dumps({'frozen_ids':len(ids),'preexisting_valid':len(valid),'missing':len(missing),'missing_ids':missing},indent=2))
unexpected=set(missing)-exp
if unexpected: raise SystemExit(f'[FAIL] unexpected missing IDs outside frozen v1.0.2 failure set: {sorted(unexpected)}')
if not missing: print('[PASS] all 61 already present; recovery fetch skipped')
(E/'raw/retrieval_v103/MISSING_IDS_v103.txt').write_text('\n'.join(missing)+('\n' if missing else ''))
PY

MISS=$(grep -c . raw/retrieval_v103/MISSING_IDS_v103.txt || true)
if [ "$MISS" -gt 0 ]; then
  python "$FIX/scripts/recover_missing_v103.py" \
    --ids raw/retrieval_v103/MISSING_IDS_v103.txt \
    --outdir raw/genbank \
    --manifest raw/retrieval_v103/STAGE05E3A_RECOVERY_MANIFEST_v103.tsv \
    --diagnostics raw/retrieval_v103/STAGE05E3A_RECOVERY_DIAGNOSTICS_v103.tsv \
    2>&1 | tee logs/recover_missing_v103.log
fi

python - "$E3A" <<'PY'
import csv, json, hashlib, sys
from pathlib import Path
from Bio import SeqIO
E=Path(sys.argv[1]); ids=[x.strip() for x in (E/'frozen/STAGE05E3A_SEQUENCE_IDS.txt').read_text().splitlines() if x.strip()]
rows=[]
v102={}
p=E/'raw/STAGE05E3A_ENA_RETRIEVAL_MANIFEST_v102.tsv'
if p.exists():
 with p.open() as h:
  for r in csv.DictReader(h,delimiter='\t'): v102[r['sequence_id']]=r
v103={}
p=E/'raw/retrieval_v103/STAGE05E3A_RECOVERY_MANIFEST_v103.tsv'
if p.exists():
 with p.open() as h:
  for r in csv.DictReader(h,delimiter='\t'): v103[r['sequence_id']]=r
for sid in ids:
 f=E/'raw/genbank'/f'{sid}.gbff'
 if not f.exists(): raise SystemExit(f'[FAIL] missing after recovery: {sid}')
 rs=list(SeqIO.parse(str(f),'genbank'))
 if len(rs)!=1 or len(rs[0].seq)<=0: raise SystemExit(f'[FAIL] invalid GenBank after recovery: {sid}')
 b=f.read_bytes(); src=v103.get(sid) or v102.get(sid) or {}
 rows.append({'sequence_id':sid,'status':'VALID','bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),
              'provider':src.get('provider','PREEXISTING_VALIDATED_GBFF'),'selected_route':src.get('selected_route','PREEXISTING_VALIDATED_GBFF'),
              'source_format':src.get('source_format','GenBank'),'output_format':src.get('output_format','GenBank')})
out=E/'raw/STAGE05E3A_COMBINED_PROVIDER_MANIFEST_v103.tsv'
with out.open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
prov={'stage':'Stage05E3A','transport_version':'v1.0.3','n_frozen_ids':len(ids),'n_validated_genbank_records':len(rows),
      'provider_strategy':'Preserve v1.0.2 ENA-valid records; recover RefSeq NC_* through EMBL-EBI Dbfetch refseqn; retry remaining INSDC through ENA Dbfetch.',
      'scientific_universe_changed':False,'family_membership_changed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,
      'downstream_analysis_script_changed':False}
(E/'outputs/STAGE05E3A_RETRIEVAL_PROVENANCE_v103.json').write_text(json.dumps(prov,indent=2)+'\n')
print(json.dumps(prov,indent=2))
PY

echo '[PASS] 61/61 GenBank-compatible records validated across hybrid EBI providers'
python scripts/analyze_official_context.py \
  --frozen frozen \
  --gbdir raw/genbank \
  --outdir outputs \
  2>&1 | tee logs/analyze_official_context_v103.log

python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/GATE_STAGE05E3A.json')
if not p.exists(): raise SystemExit('[FAIL] missing E3A gate')
x=json.loads(p.read_text())
if not x.get('gate_pass'): raise SystemExit(f'[FAIL] E3A gate: {x}')
if x.get('n_genbank_records_parsed') != 61: raise SystemExit(f"[FAIL] parsed {x.get('n_genbank_records_parsed')} records, expected 61")
print('[PASS] Stage05E3A official-context analysis complete after hybrid EBI v1.0.3 recovery')
PY
