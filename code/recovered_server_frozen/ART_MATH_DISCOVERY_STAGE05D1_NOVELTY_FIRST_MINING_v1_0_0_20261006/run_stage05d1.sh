#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
THREADS="${1:-8}"
REF="$ROOT/references/broad_rt"
WORK="$ROOT/work"; OUT="$ROOT/outputs"; LOG="$ROOT/logs"
mkdir -p "$WORK" "$OUT" "$LOG"
command -v hmmscan >/dev/null || { echo '[FAIL] hmmscan missing'; exit 2; }
command -v hmmpress >/dev/null || { echo '[FAIL] hmmpress missing'; exit 2; }
command -v blastp >/dev/null || { echo '[FAIL] blastp missing'; exit 2; }
command -v blastn >/dev/null || { echo '[FAIL] blastn missing'; exit 2; }
command -v makeblastdb >/dev/null || { echo '[FAIL] makeblastdb missing'; exit 2; }
[ -s "$REF/RVT-All.hmm" ] && [ -s "$REF/cdd-pfamA.hmm" ] || { echo '[FAIL] run ./scripts/prepare_references.sh first'; exit 3; }
# Freeze reference bytes for this run.
sha256sum -c "$REF/REFERENCE_SHA256.txt"
# HMM databases
[ -f "$REF/RVT-All.hmm.h3m" ] || hmmpress "$REF/RVT-All.hmm" >/dev/null
[ -f "$REF/cdd-pfamA.hmm.h3m" ] || hmmpress "$REF/cdd-pfamA.hmm" >/dev/null
# myRT primary focused scan, using the published explicit thresholds (-E/--domE 1e-4)
hmmscan --cpu "$THREADS" -E 1e-4 --domE 1e-4 --noali \
  --tblout "$WORK/myrt.tbl" --domtblout "$WORK/myrt.domtbl" \
  "$REF/RVT-All.hmm" "$ROOT/inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa" > "$LOG/myrt.log"
# Extract primary-hit proteins for refinement
python - "$WORK/myrt.domtbl" "$ROOT/inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa" "$WORK/myrt_primary_hits.faa" <<'PY'
import sys
from Bio import SeqIO
hits=set()
for l in open(sys.argv[1]):
 if l.startswith('#') or not l.strip(): continue
 p=l.split()
 if float(p[6])<=1e-4 and float(p[12])<=1e-4: hits.add(p[3])
with open(sys.argv[3],'w') as o:
 for r in SeqIO.parse(sys.argv[2],'fasta'):
  if r.id in hits: SeqIO.write(r,o,'fasta')
print('primary_hit_unique_targets',len(hits))
PY
# Expanded CDD/Pfam refinement only for primary hits.
if [ -s "$WORK/myrt_primary_hits.faa" ]; then
 hmmscan --cpu "$THREADS" -E 1e-3 --domE 1e-3 --noali \
   --tblout "$WORK/cdd.tbl" --domtblout "$WORK/cdd.domtbl" \
   "$REF/cdd-pfamA.hmm" "$WORK/myrt_primary_hits.faa" > "$LOG/cdd.log"
else
 : > "$WORK/cdd.tbl"; : > "$WORK/cdd.domtbl"
fi
# Registered-reference proximity is descriptive only, not complete Anthropic catalogue dereplication.
makeblastdb -in "$ROOT/references/REGISTERED_ART_RT_REFERENCES.faa" -dbtype prot -out "$WORK/art_rt_db" >/dev/null
makeblastdb -in "$ROOT/references/REGISTERED_PARTNER_REFERENCES.faa" -dbtype prot -out "$WORK/partner_db" >/dev/null
makeblastdb -in "$ROOT/references/REGISTERED_ARRAY_CORE_REFERENCES.fna" -dbtype nucl -out "$WORK/array_db" >/dev/null
blastp -query "$ROOT/inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa" -db "$WORK/art_rt_db" -num_threads "$THREADS" -evalue 1e-5 \
 -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' -out "$WORK/rt_vs_registered.tsv"
blastp -query "$ROOT/inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa" -db "$WORK/partner_db" -num_threads "$THREADS" -evalue 1e-5 \
 -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' -out "$WORK/allproteins_vs_partner.tsv"
blastn -task blastn-short -query "$ROOT/inputs/NOVELTY_FIRST_REPEAT_CORES.fna" -db "$WORK/array_db" -num_threads "$THREADS" -evalue 1e-3 \
 -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' -out "$WORK/cores_vs_registered.tsv"
python "$ROOT/scripts/analyze_novelty.py" --root "$ROOT"
