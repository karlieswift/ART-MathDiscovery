#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"; T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
mkdir -p "$T/e3c/outputs" "$T/e3a/frozen" "$T/e3a/outputs" "$T/out" "$T/work"
cat > "$T/e3c/outputs/STAGE05E3C_REVIEW_QUEUE.tsv" <<'TSV'
family_id	e3c_review_tier	official_supported_witness_n	official_witness_median_aa
F1	E3C1_PRIORITY_BIOLOGICAL_CLOSURE	1	5
F2	E3C1_PRIORITY_BIOLOGICAL_CLOSURE	1	6
F3	E3C1_PRIORITY_BIOLOGICAL_CLOSURE	1	7
F4	E3C5_KNOWN_SYSTEM_OR_SPECIAL_CONFLICT_REVIEW	0	
TSV
cat > "$T/e3a/frozen/STAGE05E3A_E2A_WITNESS_PAIRS.tsv" <<'TSV'
family_id	candidate_key	sequence_id	repeat_group_id	repeat_assay	signature_repeat_core	signature_FAM50_cluster	record_key	aa_sha256
F1	A1	S1	R1	exact	ACGTACGT	C1	K1	h1
F2	A2	S2	R2	exact	AACCGGTT	C2	K2	h2
F3	A3	S3	R3	exact	AAAACCCC	C3	K3	h3
F4	A4	S4	R4	exact	CCCCGGGG	C4	K4	h4
TSV
cat > "$T/e3a/outputs/STAGE05E3A_CANDIDATE_OFFICIAL_CONTEXT.tsv" <<'TSV'
family_id	candidate_key	sequence_id	has_exact_primary_joint_witness	witness_official_match_class	witness_official_protein_id	witness_official_product
F1	A1	S1	True	HIGH_OVERLAP	P1.1	hypothetical
F2	A2	S2	True	HIGH_OVERLAP	P2.1	hypothetical
F3	A3	S3	True	NEAR_EXACT	P3.1	hypothetical
F4	A4	S4	True			
TSV
cat > "$T/e3a/outputs/STAGE05E3A_NEIGHBORHOOD_OFFICIAL_PROTEINS.faa" <<'FA'
>P1.1
AAAAA
>P2.1
CCCCCC
>P3.1
DDDDDDD
>NEIGHBOR_NOT_WITNESS
MMMMMMMMMMMMMMMMMMMM
FA
python "$HERE/scripts/stage05e3d.py" --e3c "$T/e3c" --e3a "$T/e3a" --project-root "$T" --outdir "$T/out" --workdir "$T/work" --expected-families 4 >/dev/null
python - <<PY
import json,csv
x=json.load(open('$T/out/STAGE05E3D_SUMMARY.json')); assert x['gate_pass']; assert x['total_official_witness_materialized_n']==3; assert x['neighborhood_proteins_materialized'] is False
r=list(csv.DictReader(open('$T/out/STAGE05E3D_FAMILY_SEQUENCE_SUMMARY.tsv'),delimiter='\t'))
assert [int(float(z['official_witness_record_n'])) for z in r]==[1,1,1,0]
assert 'NEIGHBOR_NOT_WITNESS' not in open('$T/out/family_sequences/F1.official_witness.faa').read()
PY
echo '[PASS] Stage05E3D v1.0.2 exact-witness materialization selftest'
