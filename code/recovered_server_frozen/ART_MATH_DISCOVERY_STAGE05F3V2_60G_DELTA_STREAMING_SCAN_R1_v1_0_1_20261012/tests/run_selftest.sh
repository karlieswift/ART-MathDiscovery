#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd); T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
python -m py_compile "$S/scripts/"*.py
for x in "$S"/*.sh "$S/scripts/"*.sh; do bash -n "$x"; done
(cd "$S" && sha256sum -c tests/R1_UNCHANGED_SCIENCE_SHA256.txt) >/dev/null

# Original translation invariant.
cat > "$T/x.fa" <<'EOF'
>ctg1
ATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAATAA
EOF
python "$S/scripts/translate_sixframe.py" --accession TEST --min-aa 45 --stats "$T/stats.json" < "$T/x.fa" > "$T/o.faa"
grep -q '^>TEST|ctg1|+|f0|' "$T/o.faa"

# R1 full-order selector must reproduce the exact v1.0.0 top-K prefix.
cat > "$T/meta.csv" <<'EOF'
accession,assay_type,organism,tax_id,scientific_name
SRR100001,WGS,viral metagenome,1,viral metagenome
SRR100002,WGS,phage virome,2,phage virome
SRR100003,WGS,virus sample,3,virus sample
SRR100004,METAGENOMIC,marine metagenome,4,marine metagenome
SRR100005,METAGENOMIC,soil metagenome,5,soil metagenome
SRR100006,METAGENOMIC,human gut metagenome,6,human gut metagenome
SRR100007,METAGENOMIC,mouse intestinal metagenome,7,mouse intestinal metagenome
SRR100008,METAGENOMIC,hydrothermal metagenome,8,hydrothermal metagenome
SRR100009,METAGENOMIC,hypersaline metagenome,9,hypersaline metagenome
SRR100010,RNA-Seq,viral metagenome,10,viral metagenome
ERR200001,WGS,bacteriophage,11,bacteriophage
DRR300001,METAGENOMIC,wastewater metagenome,12,wastewater metagenome
EOF
: > "$T/exclude.txt"
python "$S/scripts/select_logan_metadata_v100_reference.py" --exclude "$T/exclude.txt" --k 2 --out "$T/v100.tsv" < "$T/meta.csv" 2>/dev/null
mkdir -p "$T/raw" "$T/new"
python "$S/scripts/stream_logan_candidates_r1.py" --exclude "$T/exclude.txt" --out-dir "$T/raw" < "$T/meta.csv" 2>/dev/null
for X in A B C D; do
  if [[ -s "$T/raw/$X.unsorted.tsv" ]]; then
    LC_ALL=C sort -t $'\t' -k2,2n -k3,3 "$T/raw/$X.unsorted.tsv" > "$T/$X.sorted.tsv"
    python "$S/scripts/assign_candidate_rank_r1.py" "$T/$X.sorted.tsv" "$T/new/F3V2_DELTA_CANDIDATES_$X.tsv"
  else
    printf 'stratum\tdeterministic_rank\taccession\tselection_hash64\tassay_type\torganism\ttax_id\tscientific_name\n' > "$T/new/F3V2_DELTA_CANDIDATES_$X.tsv"
  fi
done
python "$S/scripts/verify_candidate_prefix_r1.py" "$T/v100.tsv" "$T/new" "$T/prefix.json" >/dev/null
python - "$T/prefix.json" <<'PY'
import json,sys
j=json.load(open(sys.argv[1])); assert j['gate_pass']; assert all(v['exact_prefix'] for v in j['strata'].values())
PY

# Basic FASTA counter invariant.
printf '>a\nACGT\n>b\nAACCGG\n' | python "$S/scripts/count_fasta_stream.py" | grep -q '"bp": 10'
echo '[PASS] F3V2 R1 syntax + translation + exact v1.0.0 candidate-prefix selftest'
