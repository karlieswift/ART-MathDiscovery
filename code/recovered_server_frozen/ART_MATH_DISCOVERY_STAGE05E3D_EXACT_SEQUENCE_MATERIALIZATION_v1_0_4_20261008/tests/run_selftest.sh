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
# The second record intentionally mentions P1.1 in provenance metadata.  v1.0.3
# tokenized the whole header and could therefore mark P1.1 ambiguous.  v1.0.4
# must use only the structured first token for identity.
cat > "$T/e3a/outputs/STAGE05E3A_NEIGHBORHOOD_OFFICIAL_PROTEINS.faa" <<'FA'
>FASTA:aaa|S1|P1.1|hypothetical witness=P9.9
AAAAA
>FASTA:bbb|S9|P9.9|other witness=P1.1
QQQQQ
>FASTA:ccc|S2|P2.1|hypothetical
CCCCCC
>FASTA:ddd|S3|P3.1|hypothetical
DDDDDDD
FA
python "$HERE/scripts/stage05e3d.py" --e3c "$T/e3c" --e3a "$T/e3a" --project-root "$T" --outdir "$T/out" --workdir "$T/work" --expected-families 4 >/dev/null
python - <<PY
import json,csv
x=json.load(open('$T/out/STAGE05E3D_SUMMARY.json')); assert x['gate_pass']; assert x['total_official_witness_materialized_n']==3
r=list(csv.DictReader(open('$T/out/STAGE05E3D_OFFICIAL_WITNESS_AUDIT.tsv'),delimiter='\t'))
assert [z['witness_official_protein_id'] for z in r]==['P1.1','P2.1','P3.1']
assert r[0]['resolution_route']=='sequence_id+protein_id'
assert open('$T/out/family_sequences/F1.official_witness.faa').read().count('QQQQQ')==0
PY
echo '[PASS] Stage05E3D v1.0.4 structured-header/candidate-specific witness lookup selftest'
