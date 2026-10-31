#!/usr/bin/env bash
set -euo pipefail
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
mkdir -p "$T/f1/outputs" "$T/run"
cat > "$T/f1/outputs/STAGE05F1_RESCUE_PRIORITY.tsv" <<'EOF'
query_id	aa_sha256	aa_length	genome_n	candidate_n	instance_n	prodigal_relations	yxdd_like	formal_rvt_profile	exploratory_rvt_profile	best_rvt_target	best_full_evalue	best_domain_i_evalue	rescue_class	candidate_keys	accessions
SF_1	h1	100	1	1	1	NO_SAME_STRAND_PRODIGAL_OVERLAP	True	True	True	RVT-X	1e-8	1e-8	FORMAL_PROFILE_GENE_CALLER_RESCUE	A::R1	A
SF_2	h2	100	2	2	2	CONTAINED_IN_EXISTING_PRODIGAL_ORF;NO_SAME_STRAND_PRODIGAL_OVERLAP	False	True	True	RVT-Y	1e-7	1e-7	FORMAL_PROFILE_GENE_CALLER_RESCUE	B::R1;C::R1	B;C
EOF
cat > "$T/f1/outputs/STAGE05F1_SIXFRAME_PEPTIDE_INSTANCES.tsv" <<'EOF'
peptide_instance_id	sequence_id	start0	end0	strand	frame	aa_length	aa_sha256	edge_partial	prodigal_relation	same_strand_prodigal_overlap_n	candidate_keys
p1	A	100	400	+	0	100	h1	False	NO_SAME_STRAND_PRODIGAL_OVERLAP	0	A::R1
p2	B	100	400	+	0	100	h2	False	CONTAINED_IN_EXISTING_PRODIGAL_ORF	1	B::R1
p3	C	100	400	+	0	100	h2	False	NO_SAME_STRAND_PRODIGAL_OVERLAP	0	C::R1
EOF
cat > "$T/f1/outputs/STAGE05F1_RVT_ALL.domtblout" <<'EOF'
#
RVT-X - 200 SF_1 - 100 1e-8 50 0 1 1 1e-8 1e-8 45 0 20 150 10 90 8 95 0.95 desc
RVT-Y - 250 SF_2 - 100 1e-7 40 0 1 1 1e-7 1e-7 35 0 30 170 5 80 3 85 0.90 desc
EOF
echo '{"gate":"PASS_STAGE05F1_ANNOTATION_INDEPENDENT_SIXFRAME_RT_RESCUE","gate_pass":true}' > "$T/f1/outputs/STAGE05F1_GATE.json"
cd "$T/run"
python "$OLDPWD/scripts/stage05f1r1.py" "$T/f1" . --outdir outputs >/dev/null
python - <<'PYCHECK'
import json
x=json.load(open('outputs/STAGE05F1R1_GATE.json'))
assert x['f1_formal_gene_caller_rescue_unique_aa_n']==2
assert x['strict_formal_gene_model_rescue_unique_aa_n']==1
assert x['mixed_fanout_formal_unique_aa_n']==1
print('[PASS] Stage05F1R1 strict fanout/instance selftest')
PYCHECK
