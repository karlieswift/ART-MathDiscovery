#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
mkdir -p "$T/inp/inputs" "$T/s5"
cat > "$T/inp/inputs/ART_EXACT_SOURCE_FASTA_PROVENANCE.tsv" <<'TSV'
panel_id	author_label	accession	source_mode	source_origin	fasta_path	header	length_bp	sha256	acgt_only	single_record	truth_coordinates_used	repeat_motif_used
TSV
for acc in OR770614.1 OR836606.1 LC680885.1 MN091626.1 MZ779063.1; do
  printf '>%s\nACGTACGTACGT\n' "$acc" > "$T/inp/inputs/$acc.frozen.fasta"
  h=$(sha256sum "$T/inp/inputs/$acc.frozen.fasta" | awk '{print $1}')
  printf 'P\tX\t%s\tLOCAL\tX\tX\t%s\t12\t%s\tTrue\tTrue\tFalse\tFalse\n' "$acc" "$acc" "$h" >> "$T/inp/inputs/ART_EXACT_SOURCE_FASTA_PROVENANCE.tsv"
  mkdir -p "$T/s5/outputs/genomes/$acc"
  cat > "$T/s5/outputs/genomes/$acc/STAGE03A_UNION_REGIONS.tsv" <<EOF
sequence_id	region_id	region_rank	support_families
$acc	R1	1	ACNV_FAMILY;GENERIC_RECURRENCE
$acc	R2	2	ACNV_FAMILY
EOF
  cat > "$T/s5/outputs/genomes/$acc/STAGE03A_SELECTED_PROPOSAL_WINDOWS.tsv" <<EOF
sequence_id	proposal_id	family
$acc	A1	ACNV_FAMILY
$acc	Q1	GENERIC_RECURRENCE
EOF
  cat > "$T/s5/outputs/genomes/$acc/COMPLETION.json" <<EOF
{"sequence_id":"$acc","complete":true}
EOF
done
python "$HERE/scripts/recover_stage05b_panel.py" --project-root "$T" --exact-inputs "$T/inp" --stage05b "$T/s5" --outdir "$T/out" --workdir "$T/work" >/dev/null
python - <<PY
import json,csv
x=json.load(open('$T/out/ART_PANEL_STAGEA1_GATE.json')); assert x['gate_pass'] and x['n_recovered']==5
r=list(csv.DictReader(open('$T/out/ART_PANEL_STAGE05B_SCORE_RECOVERY_SUMMARY.tsv'),delimiter='\t'))
assert len(r)==5 and all(int(x['dual_arm_regions_n'])==1 for x in r)
PY
echo '[PASS] ART panel existing-score recovery selftest'
