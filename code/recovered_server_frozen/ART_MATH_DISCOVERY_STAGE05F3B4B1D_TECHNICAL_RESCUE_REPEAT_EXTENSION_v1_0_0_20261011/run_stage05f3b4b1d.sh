#!/usr/bin/env bash
set -euo pipefail
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE"/{outputs,logs,frozen}
ER1="$P/ART_MATH_DISCOVERY_STAGE05F3B3ER1_FOLDSEEK_ONLY_TECHNICAL_REPAIR_v1_0_0_20261011/outputs/F3B3ER1_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv"
R4="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_v1_0_0_20261011/outputs"
FM="$R4/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv"
IM="$R4/F3B1R3_INSTANCE_MAP_WITH_CANONICAL_SOURCE.tsv"
S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010/work"
B1="$P/ART_MATH_DISCOVERY_STAGE05F3B4B1_FAMILY_WIDE_REPEAT_ARRAY_ARCHITECTURE_v1_0_0_20261011"
B1C="$P/ART_MATH_DISCOVERY_STAGE05F3B4B1C_SHORTLIST_WEAK_REPEAT_CONTROL_v1_0_0_20261011"
SV="$B1/outputs/F3B4B1_FAMILY_REPEAT_ARRAY_SUMMARY_16.tsv"
WK="$B1C/outputs/F3B4B1C_WEAK_FAMILY_REPEAT_SUMMARY_9.tsv"
UPDET="$B1/scripts/repeat_array_scan.py"
UPPLC="$P/ART_MATH_DISCOVERY_STAGE05F3B4B0_FAMILY_WIDE_CONTEXT_COMPLETENESS_RESCUE_v1_0_0_20261011/scripts/family_wide_context_rescue.py"
for f in "$ER1" "$FM" "$IM" "$S0/F3B1_orfs_aa" "$S0/F3B1_orfs_aa.index" "$S0/genomeDB" "$S0/genomeDB.index" "$SV" "$WK" "$UPDET" "$UPPLC"; do [[ -e "$f" ]] || { echo "[FATAL] missing $f"; exit 2; }; done
"$HERE/tests/run_selftest.sh"

echo '[1/3] freeze all-five repaired-family architecture extension before repeat results'
python - "$HERE" "$ER1" "$FM" "$IM" "$SV" "$WK" "$UPDET" "$UPPLC" <<'PY'
import csv,json,hashlib,sys
from pathlib import Path
here,er1,fm,im,sv,wk,updet,upplc=map(Path,sys.argv[1:])
def sha(p):
 h=hashlib.sha256();h.update(p.read_bytes());return h.hexdigest()
r=list(csv.DictReader(open(er1),delimiter='\t'))
t=sorted([x for x in r if x.get('b3er1_resolution_source')=='FOLDSEEK_ONLY_TECHNICAL_REPAIR'],key=lambda x:int(x['shortlist_rank']))
assert len(r)==32 and len(t)==5
assert sum(x['b3er1_final_best_tier'] in {'STRONG','SUGGESTIVE'} for x in t)==2
assert sum(x['b3er1_final_best_tier']=='WEAK' for x in t)==3
# vendored detector/placement helpers must be byte-identical to upstream frozen implementations
assert sha(here/'scripts/frozen_b4b1_repeat_array_scan.py')==sha(updet)
assert sha(here/'scripts/frozen_b4b0_family_wide_context_rescue.py')==sha(upplc)
lock={'stage':'Stage05F3B4B1D-implementation-freeze','gate':'PASS_F3B4B1D_TECHNICAL_RESCUE_ARCHITECTURE_EXTENSION_FROZEN_PRE_RESULTS','gate_pass':True,
 'analysis_role':'POST_HOC_TECHNICAL_RESCUE_EXTENSION_ALL_AND_ONLY_FIVE_REPAIRED_FAMILIES',
 'target_rule':'all and only five F3B3E-R1 families with b3er1_resolution_source=FOLDSEEK_ONLY_TECHNICAL_REPAIR; no selection by repeat result',
 'target_family_n':5,'target_final_strong_or_suggestive_family_n':2,'target_final_weak_family_n':3,
 'target_family_ids':[x['family_rep_id'] for x in t],'target_candidate_ids':[x['candidate_id'] for x in t],
 'repeat_detector':'byte-identical vendored F3B4B1 detector: k=8,10,12,14,16; >=3 copies; start-gap 60-360; exact-core <=80; same strict/broad annotations',
 'placement_engine':'byte-identical vendored F3B4B0 six-frame exact placement helper',
 'original_b4b1_16_membership_changed':False,'original_b4b1c_weak9_control_changed':False,'two_resource_limited_unresolved_treated_as_negative':False,
 'threshold_retuning':False,'novelty_used':False,'partner_analysis_performed':False,
 'er1_summary_sha256':sha(er1),'r4_family_map_sha256':sha(fm),'r4_instance_map_sha256':sha(im),'original_b4b1_summary_sha256':sha(sv),'original_b4b1c_summary_sha256':sha(wk),
 'upstream_b4b1_detector_sha256':sha(updet),'vendored_b4b1_detector_sha256':sha(here/'scripts/frozen_b4b1_repeat_array_scan.py'),
 'upstream_b4b0_placement_sha256':sha(upplc),'vendored_b4b0_placement_sha256':sha(here/'scripts/frozen_b4b0_family_wide_context_rescue.py'),
 'extension_script_sha256':sha(here/'scripts/run_extension.py'),
 'claim_boundary':'Architecture analysis of the five newly structurally resolved families only. Integrated 18-vs-12 comparison is descriptive/post-discovery; original frozen analyses remain preserved.'}
(here/'frozen/F3B4B1D_IMPLEMENTATION_LOCK.json').write_text(json.dumps(lock,indent=2)+'\n')
print(json.dumps(lock,indent=2))
PY

echo '[2/3] exact-place all family members and scan repeats with unchanged frozen detector'
python "$HERE/scripts/run_extension.py" \
  "$ER1" "$FM" "$IM" "$S0/F3B1_orfs_aa" "$S0/F3B1_orfs_aa.index" "$S0/genomeDB" "$S0/genomeDB.index" \
  "$SV" "$WK" "$HERE/outputs" "$HERE/frozen/F3B4B1D_IMPLEMENTATION_LOCK.json" \
  "$HERE/scripts/frozen_b4b1_repeat_array_scan.py" "$HERE/scripts/frozen_b4b0_family_wide_context_rescue.py" 5 \
  2>&1 | tee "$HERE/logs/f3b4b1d_extension.log"

echo '[3/3] hashes'
(cd "$HERE" && sha256sum frozen/F3B4B1D_IMPLEMENTATION_LOCK.json outputs/* > outputs/STAGE05F3B4B1D_SHA256.tsv)
echo PASS_STAGE05F3B4B1D_TECHNICAL_RESCUE_REPEAT_EXTENSION_RUN_COMPLETE
