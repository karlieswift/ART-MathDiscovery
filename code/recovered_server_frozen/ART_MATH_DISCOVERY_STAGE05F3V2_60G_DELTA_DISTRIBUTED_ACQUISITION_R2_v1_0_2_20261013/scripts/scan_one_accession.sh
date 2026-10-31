#!/usr/bin/env bash
set -euo pipefail
S=$1; ACC=$2; HMM=$3; CPU=$4
Z="$S/data/contigs/$ACC.contigs.fa.zst"; D="$S/work/per_accession/$ACC"; mkdir -p "$D" "$S/outputs/formal_faa" "$S/logs/per_accession" "$S/work/done"
[[ -s "$Z" ]] || { echo "missing $Z" >&2; exit 2; }
[[ -e "$S/work/done/$ACC.done" ]] && exit 0
FAA="$D/$ACC.sixframe.faa"
zstd -dc "$Z" | python "$S/scripts/translate_sixframe.py" --accession "$ACC" --min-aa 45 --stats "$D/translation_stats.json" > "$FAA"
hmmscan --cpu "$CPU" -E 0.01 --domE 0.01 --tblout "$D/RVT_ALL.tblout" --domtblout "$D/RVT_ALL.domtblout" "$HMM" "$FAA" > "$D/RVT_ALL.stdout" 2> "$D/RVT_ALL.stderr"
python "$S/scripts/parse_hmmscan.py" "$D/RVT_ALL.domtblout" "$D/candidates.tsv"
python "$S/scripts/extract_fasta_ids.py" "$D/candidates.tsv" "$FAA" "$S/outputs/formal_faa/$ACC.formal.faa"
if [[ "${KEEP_ORF_FASTA:-0}" != 1 ]]; then rm -f "$FAA"; fi
touch "$S/work/done/$ACC.done"
echo "[PASS] $ACC"
