#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
R4="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_v1_0_0_20261011"
A="$P/ART_MATH_DISCOVERY_STAGE05F3B3A_POSTFREEZE_PRIORITY_AND_CATALYTIC_TOPOLOGY_v1_0_0_20261011"
B="$P/ART_MATH_DISCOVERY_STAGE05F3B3B_KNOWN_SYSTEM_AUDIT_AND_STRUCTURE_HANDOFF_v1_0_0_20261011"
D="$P/ART_MATH_DISCOVERY_STAGE05F3B3D_EXACT_PROTOCOL_TECHNICAL_RETRY_v1_0_0_20261011"
S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010"
mkdir -p "$S"/{frozen,outputs,logs}

R4MAP="$R4/outputs/F3B1R3_INSTANCE_MAP_WITH_CANONICAL_SOURCE.tsv"
DSUM="$D/outputs/F3B3D_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv"
DGATE="$D/outputs/STAGE05F3B3D_GATE.json"
QFA="$A/outputs/F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.faa"
BAUD="$B/outputs/F3B3B_KNOWN_SYSTEM_AUDIT_32.tsv"
GDB="$S0/work/genomeDB"
GIDX="$S0/work/genomeDB.index"
LOCK="$S/frozen/F3B4A_IMPLEMENTATION_LOCK.json"

for x in python sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for x in "$R4MAP" "$DSUM" "$DGATE" "$QFA" "$BAUD" "$GDB" "$GIDX" "$S/scripts/reconstruct_survivor_loci.py"; do [[ -s "$x" ]] || { echo "[FATAL] missing input: $x" >&2; exit 2; }; done

echo '[0/4] selftest'
(cd "$S" && ./tests/run_selftest.sh)

if [[ ! -s "$LOCK" ]]; then
  echo '[1/4] freeze exact-locus reconstruction implementation BEFORE locus results'
  python - "$R4MAP" "$DSUM" "$DGATE" "$QFA" "$BAUD" "$S/scripts/reconstruct_survivor_loci.py" "$S/run_stage05f3b4a.sh" "$LOCK" <<'PY'
import hashlib,json,sys,csv
r4,ds,dg,qfa,ba,py,run,out=sys.argv[1:]
def h(p):
 x=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(16<<20),b''):x.update(b)
 return x.hexdigest()
g=json.load(open(dg));assert g['gate_pass'] is True
with open(ds) as f:r=list(csv.DictReader(f,delimiter='\t'))
s=[x for x in r if x['final_best_tier'] in ('STRONG','SUGGESTIVE')]
assert len(r)==32 and len(s)==16
lock={
 'stage':'Stage05F3B4A-implementation-freeze','gate':'PASS_F3B4A_NATIVE_LOCUS_RECONSTRUCTION_FROZEN_PRE_LOCUS','gate_pass':True,
 'survivor_rule':'all and only F3B3D rows with final_best_tier in {STRONG,SUGGESTIVE}','survivor_expected_n':16,
 'flank_bp_each_side':10000,'translation_table':11,'start_codon_requirement':'none',
 'placement_rule':'unique exact peptide occurrence in six-frame translation of canonical source contig; 0 or >1 exact placements cause HOLD; no manual tie-break',
 'candidate_orientation_output':True,'genomic_orientation_output':True,
 'annotation_coordinates_used':False,'repeat_or_partner_analysis_performed_in_this_stage':False,'novelty_decision_performed':False,
 'f3b3d_summary_sha256':h(ds),'f3b3d_gate_sha256':h(dg),'r4_instance_map_sha256':h(r4),'shortlist_fasta_sha256':h(qfa),'b3b_audit_sha256':h(ba),
 'implementation_sha256':h(py),'runner_sha256':h(run),
 'claim_boundary':'Exact source-locus/context reconstruction only; no ART or novelty claim.'}
with open(out,'w') as f:json.dump(lock,f,indent=2);f.write('\n')
print(json.dumps(lock,indent=2))
PY
  sha256sum "$LOCK" "$S/scripts/reconstruct_survivor_loci.py" "$S/run_stage05f3b4a.sh" > "$S/frozen/F3B4A_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$S/scripts/reconstruct_survivor_loci.py" "$S/run_stage05f3b4a.sh"
else
  echo '[1/4] [SKIP] implementation lock exists'
fi

echo '[2/4] exact native-locus reconstruction for frozen 16 structural survivors'
python "$S/scripts/reconstruct_survivor_loci.py" \
  "$R4MAP" "$DSUM" "$DGATE" "$QFA" "$BAUD" "$GDB" "$GIDX" \
  "$S/outputs" "$LOCK" 10000 \
  2>&1 | tee "$S/logs/f3b4a_reconstruct.log"

echo '[3/4] hashes'
find "$S/outputs" -maxdepth 2 -type f -print0 | sort -z | xargs -0 sha256sum > "$S/outputs/STAGE05F3B4A_SHA256.tsv"

echo '[4/4] compact result'
cat "$S/outputs/STAGE05F3B4A_GATE.json"
echo
echo '=== SURVIVOR LOCUS MAP ==='
column -t -s $'\t' "$S/outputs/F3B4A_SURVIVOR_NATIVE_LOCUS_MAP.tsv" | head -17 || true
echo
echo 'PASS_STAGE05F3B4A_SURVIVOR_NATIVE_LOCUS_RECONSTRUCTION_RUN_COMPLETE'
