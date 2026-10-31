#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then echo "usage: $0 <E3A_ROOT>" >&2; exit 2; fi
E3A="$(cd "$1" && pwd)"
FIX="$(cd "$(dirname "$0")" && pwd)"
cd "$E3A"
./run_selftest.sh
"$FIX/run_selftest.sh"
mkdir -p raw/genbank raw/retrieval_v104 outputs logs

python - "$E3A" <<'PY'
import json, sys
from pathlib import Path
from Bio import SeqIO
E=Path(sys.argv[1]); ids=[x.strip() for x in (E/'frozen/STAGE05E3A_SEQUENCE_IDS.txt').read_text().splitlines() if x.strip()]
if len(ids)!=61: raise SystemExit(f'[FAIL] frozen sequence count {len(ids)} != 61')
valid=[]; missing=[]
for sid in ids:
 p=E/'raw/genbank'/f'{sid}.gbff'; ok=False
 if p.exists():
  try:
   rs=list(SeqIO.parse(str(p),'genbank')); ok=(len(rs)==1 and len(rs[0].seq)>0)
  except Exception: ok=False
 (valid if ok else missing).append(sid)
print(json.dumps({'frozen_ids':len(ids),'preexisting_valid':len(valid),'missing':len(missing),'missing_ids':missing},indent=2))
if missing not in ([], ['OK905446.1']):
 raise SystemExit(f'[FAIL] v1.0.4 expected only OK905446.1 missing; observed {missing}')
PY

if [ ! -s raw/genbank/OK905446.1.gbff ]; then
  python "$FIX/scripts/recover_ok905446_ddbj_v104.py" \
    --outdir raw/genbank \
    --manifest raw/retrieval_v104/STAGE05E3A_RECOVERY_MANIFEST_v104.tsv \
    --diagnostics raw/retrieval_v104/STAGE05E3A_RECOVERY_DIAGNOSTICS_v104.tsv \
    2>&1 | tee logs/recover_ok905446_ddbj_v104.log
fi

python - "$E3A" <<'PY'
import csv, json, hashlib, sys
from pathlib import Path
from Bio import SeqIO
E=Path(sys.argv[1]); ids=[x.strip() for x in (E/'frozen/STAGE05E3A_SEQUENCE_IDS.txt').read_text().splitlines() if x.strip()]
sources=[]
for pp in [
 E/'raw/retrieval_v104/STAGE05E3A_RECOVERY_MANIFEST_v104.tsv',
 E/'raw/retrieval_v103/STAGE05E3A_RECOVERY_MANIFEST_v103.tsv',
 E/'raw/STAGE05E3A_ENA_RETRIEVAL_MANIFEST_v102.tsv',
]:
 d={}
 if pp.exists():
  with pp.open() as h:
   for r in csv.DictReader(h,delimiter='\t'):
    if r.get('status')=='VALID': d[r['sequence_id']]=r
 sources.append(d)
rows=[]
for sid in ids:
 f=E/'raw/genbank'/f'{sid}.gbff'
 if not f.exists(): raise SystemExit(f'[FAIL] missing after v1.0.4 recovery: {sid}')
 rs=list(SeqIO.parse(str(f),'genbank'))
 if len(rs)!=1 or len(rs[0].seq)<=0: raise SystemExit(f'[FAIL] invalid GenBank after v1.0.4 recovery: {sid}')
 src={}
 for d in sources:
  if sid in d: src=d[sid]; break
 b=f.read_bytes()
 rows.append({'sequence_id':sid,'status':'VALID','bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),
              'provider':src.get('provider','PREEXISTING_VALIDATED_GBFF'),
              'selected_route':src.get('selected_route','PREEXISTING_VALIDATED_GBFF'),
              'source_format':src.get('source_format','GenBank'),
              'output_format':src.get('output_format','GenBank')})
out=E/'raw/STAGE05E3A_COMBINED_PROVIDER_MANIFEST_v104.tsv'
with out.open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
prov={'stage':'Stage05E3A','transport_version':'v1.0.4','n_frozen_ids':len(ids),'n_validated_genbank_records':len(rows),
      'provider_strategy':'Preserve validated v1.0.2 ENA and v1.0.3 EBI RefSeq records; recover only OK905446.1 via DDBJ getentry (INSDC flatfile).',
      'scientific_universe_changed':False,'family_membership_changed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,
      'downstream_analysis_script_changed':False}
(E/'outputs/STAGE05E3A_RETRIEVAL_PROVENANCE_v104.json').write_text(json.dumps(prov,indent=2)+'\n')
print(json.dumps(prov,indent=2))
PY

echo '[PASS] 61/61 GenBank-compatible records validated across ENA + EBI RefSeq + DDBJ providers'
python scripts/analyze_official_context.py \
  --frozen frozen \
  --gbdir raw/genbank \
  --outdir outputs \
  2>&1 | tee logs/analyze_official_context_v104.log

python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/GATE_STAGE05E3A.json')
if not p.exists(): raise SystemExit('[FAIL] missing E3A gate')
x=json.loads(p.read_text())
if not x.get('gate_pass'): raise SystemExit(f'[FAIL] E3A gate: {x}')
if x.get('n_genbank_records_parsed') != 61: raise SystemExit(f"[FAIL] parsed {x.get('n_genbank_records_parsed')} records, expected 61")
print('[PASS] Stage05E3A official-context analysis complete after DDBJ v1.0.4 final-record recovery')
PY
