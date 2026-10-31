#!/usr/bin/env bash
set -euo pipefail
R1="${1:?usage: run_sensitive_homology.sh <Stage05D4R1-dir> [threads]}"
THREADS="${2:-32}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/outputs"; WORK="$ROOT/work"; REFS="$ROOT/refs"; LOG="$ROOT/logs"
mkdir -p "$OUT" "$WORK" "$REFS" "$LOG"
Q="$R1/outputs/STAGE05D4R1_NR90_REPRESENTATIVES.faa"
G="$R1/outputs/GATE_STAGE05D4R1.json"
T="$R1/outputs/STAGE05D4R1_NR90_CLUSTERS.tsv"
[[ -s "$Q" && -s "$G" && -s "$T" ]] || { echo "[FAIL] missing Stage05D4-R1 inputs"; exit 2; }
python - "$G" <<'PY'
import json,sys
j=json.load(open(sys.argv[1]))
assert j['gate']=='PASS_N1N2_SEQUENCE_FAMILY_COMPRESSION_FROZEN' and j['gate_pass'] is True
assert j['n_N1N2_candidates']==602
assert j['n_input_nonformal_unique_aa']==13062
print('[PASS] frozen D4-R1 gate')
PY
command -v mmseqs >/dev/null || { echo '[FAIL] mmseqs not found'; exit 3; }
command -v wget >/dev/null || { echo '[FAIL] wget not found'; exit 3; }
DBURL='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
GZ="$REFS/uniprot_sprot.current_release.fasta.gz"
FA="$REFS/uniprot_sprot.current_release.fasta"
if [[ ! -s "$GZ" ]]; then
  wget -c --tries=20 --timeout=60 "$DBURL" -O "$GZ"
fi
if [[ ! -s "$FA" || "$GZ" -nt "$FA" ]]; then
  gzip -dc "$GZ" > "$FA.tmp"
  mv "$FA.tmp" "$FA"
fi
sha256sum "$GZ" "$FA" > "$OUT/SWISSPROT_REFERENCE_SHA256.tsv"
python - "$FA" "$OUT/SWISSPROT_REFERENCE_METADATA.json" <<'PY'
from pathlib import Path
import json,sys,hashlib,datetime
p=Path(sys.argv[1]); n=0; aas=0
with p.open() as f:
    for line in f:
        if line.startswith('>'): n+=1
        else: aas += len(line.strip())
out={'source':'UniProtKB/Swiss-Prot current_release','retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sequences':n,'total_aa':aas,'fasta_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
PY
rm -rf "$WORK/mmseqs_tmp" "$WORK/queryDB"* "$WORK/sprotDB"* "$WORK/search"*
mmseqs createdb "$Q" "$WORK/queryDB"
mmseqs createdb "$FA" "$WORK/sprotDB"
mmseqs search "$WORK/queryDB" "$WORK/sprotDB" "$WORK/search" "$WORK/mmseqs_tmp" \
  --threads "$THREADS" -s 7.5 -e 0.001 --max-seqs 25 --cov-mode 0 -c 0.0
mmseqs convertalis "$WORK/queryDB" "$WORK/sprotDB" "$WORK/search" "$OUT/STAGE05D4R2_SWISSPROT_HITS.tsv" \
  --format-output 'query,target,fident,alnlen,qstart,qend,tstart,tend,qcov,tcov,evalue,bits'
python "$ROOT/scripts/summarize_hits.py" "$Q" "$FA" "$T" "$OUT/STAGE05D4R2_SWISSPROT_HITS.tsv" "$OUT"
python "$ROOT/scripts/write_gate.py" "$R1" "$OUT"
python - "$OUT" <<'PY'
from pathlib import Path
import hashlib,sys
out=Path(sys.argv[1]); rows=[]
for p in sorted(out.iterdir()):
    if p.is_file() and p.name!='STAGE05D4R2_OUTPUT_MANIFEST.tsv':
        rows.append((p.name,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()))
with (out/'STAGE05D4R2_OUTPUT_MANIFEST.tsv').open('w') as f:
    f.write('file\tbytes\tsha256\n')
    for n,b,h in rows: f.write(f'outputs/{n}\t{b}\t{h}\n')
PY
cat "$OUT/GATE_STAGE05D4R2.json"
