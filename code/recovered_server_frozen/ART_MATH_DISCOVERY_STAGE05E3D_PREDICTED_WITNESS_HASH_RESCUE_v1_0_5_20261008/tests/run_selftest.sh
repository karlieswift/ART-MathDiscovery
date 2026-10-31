#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"; T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
mkdir -p "$T/e3c/outputs" "$T/e3a/frozen" "$T/e3a/outputs" "$T/ART_MATH_DISCOVERY_STAGE05D1_SYNTH/outputs"
cat > "$T/e3c/outputs/STAGE05E3C_REVIEW_QUEUE.tsv" <<'TSV'
family_id	e3c_review_tier	official_supported_witness_n	official_witness_median_aa
F1	E3C1_PRIORITY_BIOLOGICAL_CLOSURE	1	5
F2	E3C1_PRIORITY_BIOLOGICAL_CLOSURE	1	6
F3	E3C1_PRIORITY_BIOLOGICAL_CLOSURE	1	7
F4	E3C5_KNOWN_SYSTEM_OR_SPECIAL_CONFLICT_REVIEW	0	
TSV
python - <<PY
import hashlib
seqs={'K1':'MKKKK','K2':'MCCCCC','K3':'MDDDDDD','K4':'MEEEEE'}
with open('$T/ART_MATH_DISCOVERY_STAGE05D1_SYNTH/outputs/PREDICTED_WITNESSES.faa','w') as f:
  for k,s in seqs.items(): f.write(f'>{k}\n{s}\n')
with open('$T/e3a/frozen/STAGE05E3A_E2A_WITNESS_PAIRS.tsv','w') as f:
  f.write('family_id\tcandidate_key\tsequence_id\trepeat_group_id\trepeat_assay\tsignature_repeat_core\tsignature_FAM50_cluster\trecord_key\taa_sha256\n')
  for fam,cand,sid,rg,core,cl,k in [
    ('F1','A1','S1','R1','ACGTACGT','C1','K1'),('F2','A2','S2','R2','AACCGGTT','C2','K2'),('F3','A3','S3','R3','AAAACCCC','C3','K3'),('F4','A4','S4','R4','CCCCGGGG','C4','K4')]:
    h=hashlib.sha256(seqs[k].encode()).hexdigest()
    f.write(f'{fam}\t{cand}\t{sid}\t{rg}\texact\t{core}\t{cl}\t{k}\t{h}\n')
PY
cat > "$T/e3a/outputs/STAGE05E3A_CANDIDATE_OFFICIAL_CONTEXT.tsv" <<'TSV'
family_id	candidate_key	sequence_id	has_exact_primary_joint_witness	witness_official_match_class	witness_official_protein_id	witness_official_product
F1	A1	S1	True	HIGH_OVERLAP	P1.1	hypothetical
F2	A2	S2	True	HIGH_OVERLAP	P2.1	hypothetical
F3	A3	S3	True	NEAR_EXACT	P3.1	hypothetical
F4	A4	S4	True			
TSV
# Deliberately omit the official-protein FASTA: v1.0.5 must recover the
# family-defining predicted witness by E2A aa_sha256 instead.
python "$HERE/scripts/stage05e3d.py" --e3c "$T/e3c" --e3a "$T/e3a" --project-root "$T" --outdir "$T/out" --workdir "$T/work" --expected-families 4 >/dev/null
python - <<PY
import json,csv
x=json.load(open('$T/out/STAGE05E3D_SUMMARY.json'))
assert x['gate_pass'],x
assert x['total_predicted_witness_materialized_n']==3,x
assert x['official_protein_fasta_required'] is False
r=list(csv.DictReader(open('$T/out/STAGE05E3D_PREDICTED_WITNESS_AUDIT.tsv'),delimiter='\t'))
assert [z['record_key'] for z in r]==['K1','K2','K3'],r
assert all('STAGE05D1_SYNTH' in z['source_file'] for z in r),r
PY
echo '[PASS] Stage05E3D v1.0.5 E2A-aa_sha256 predicted-witness rescue selftest'
