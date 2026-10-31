#!/usr/bin/env bash
set -euo pipefail
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/outputs" "$HERE/logs" "$HERE/frozen"
B4A="$P/ART_MATH_DISCOVERY_STAGE05F3B4A_SURVIVOR_NATIVE_LOCUS_RECONSTRUCTION_v1_0_0_20261011/outputs/F3B4A_SURVIVOR_NATIVE_LOCUS_MAP.tsv"
R4="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_v1_0_0_20261011/outputs"
FM="$R4/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv"
IM="$R4/F3B1R3_INSTANCE_MAP_WITH_CANONICAL_SOURCE.tsv"
S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010/work"
for f in "$B4A" "$FM" "$IM" "$S0/F3B1_orfs_aa" "$S0/F3B1_orfs_aa.index" "$S0/genomeDB" "$S0/genomeDB.index"; do [[ -e "$f" ]] || { echo "[FATAL] missing $f"; exit 2; }; done

echo '[1/3] freeze family-wide technical context-rescue rule before family-context results'
python - "$HERE" "$B4A" "$FM" "$IM" <<'PY'
import csv,json,hashlib,sys
from pathlib import Path
here,b4a,fm,im=map(Path,sys.argv[1:])
def sha(p):
 h=hashlib.sha256(); h.update(p.read_bytes()); return h.hexdigest()
rows=list(csv.DictReader(open(b4a),delimiter='\t'))
fams=[]
for r in sorted(rows,key=lambda x:int(x['shortlist_rank'])):
 if r['family_rep_id'] not in fams: fams.append(r['family_rep_id'])
assert len(fams)==16
lock={
 'stage':'Stage05F3B4B0-implementation-freeze','gate':'PASS_F3B4B0_FAMILY_WIDE_CONTEXT_RESCUE_RULE_FROZEN_PRE_RESULTS','gate_pass':True,
 'survivor_family_n':16,'survivor_family_ids':fams,
 'member_rule':'all and only formal F3B2-R4 candidate instances belonging to the 16 frozen F3B3D STRONG/SUGGESTIVE survivor families',
 'translation_table':11,'start_codon_requirement':'none','placement_rule':'six-frame exact peptide replay on canonical parent contig; unique exact only eligible for technical completeness',
 'context_completeness_bins_bp_each_side':[500,2000,5000,10000],
 'technical_best_context_rule':'maximize min(left,right), then total flank, then parent contig length, then candidate_id lexical among unique exact placements',
 'repeat_annotation_or_novelty_used_for_context_selection':False,'repeat_analysis_performed':False,'partner_analysis_performed':False,
 'b4a_locus_map_sha256':sha(b4a),'r4_family_map_sha256':sha(fm),'r4_instance_map_sha256':sha(im),
 'implementation_sha256':sha(here/'scripts/family_wide_context_rescue.py'),
 'claim_boundary':'Pre-architecture technical context recovery across all members of frozen survivor families; no biological filtering or ART claim.'}
(here/'frozen/F3B4B0_IMPLEMENTATION_LOCK.json').write_text(json.dumps(lock,indent=2)+'\n')
print(json.dumps(lock,indent=2))
PY

echo '[2/3] reconstruct all formal member loci and quantify family-wide context completeness'
python "$HERE/scripts/family_wide_context_rescue.py" \
  "$B4A" "$FM" "$IM" \
  "$S0/F3B1_orfs_aa" "$S0/F3B1_orfs_aa.index" \
  "$S0/genomeDB" "$S0/genomeDB.index" \
  "$HERE/outputs" "$HERE/frozen/F3B4B0_IMPLEMENTATION_LOCK.json" \
  2>&1 | tee "$HERE/logs/f3b4b0_family_context.log"

echo '[3/3] hashes'
(cd "$HERE" && sha256sum frozen/F3B4B0_IMPLEMENTATION_LOCK.json outputs/* 2>/dev/null > outputs/STAGE05F3B4B0_SHA256.tsv || true)
echo PASS_STAGE05F3B4B0_RUN_COMPLETE
