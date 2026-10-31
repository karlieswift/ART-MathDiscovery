#!/usr/bin/env bash
set -euo pipefail
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/outputs" "$HERE/logs" "$HERE/frozen"
D="$P/ART_MATH_DISCOVERY_STAGE05F3B3D_EXACT_PROTOCOL_TECHNICAL_RETRY_v1_0_0_20261011/outputs/F3B3D_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv"
R4="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_v1_0_0_20261011/outputs"
FM="$R4/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv"
IM="$R4/F3B1R3_INSTANCE_MAP_WITH_CANONICAL_SOURCE.tsv"
S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010/work"
B1="$P/ART_MATH_DISCOVERY_STAGE05F3B4B1_FAMILY_WIDE_REPEAT_ARRAY_ARCHITECTURE_v1_0_0_20261011/outputs"
SV="$B1/F3B4B1_FAMILY_REPEAT_ARRAY_SUMMARY_16.tsv"
for f in "$D" "$FM" "$IM" "$S0/F3B1_orfs_aa" "$S0/F3B1_orfs_aa.index" "$S0/genomeDB" "$S0/genomeDB.index" "$SV"; do [[ -e "$f" ]] || { echo "[FATAL] missing $f"; exit 2; }; done
"$HERE/tests/run_selftest.sh"

echo '[1/2] freeze same-shortlist WEAK comparator before control results'
python - "$HERE" "$D" "$FM" "$IM" "$SV" <<'PY'
import csv,json,hashlib,sys
from pathlib import Path
here,d,fm,im,sv=map(Path,sys.argv[1:])
def sha(p):
 h=hashlib.sha256();h.update(p.read_bytes());return h.hexdigest()
r=list(csv.DictReader(open(d),delimiter='\t'))
weak=[x for x in r if x['final_best_tier']=='WEAK'];unres=[x for x in r if x['final_best_tier'].startswith('BOLTZ_TECHNICAL_FAILURE')]
assert len(r)==32 and len(weak)==9 and len(unres)==7
lock={'stage':'Stage05F3B4B1C-implementation-freeze','gate':'PASS_F3B4B1C_WEAK_COMPARATOR_FROZEN_PRE_CONTROL_RESULTS','gate_pass':True,
 'analysis_role':'POST_DISCOVERY_SAME_SHORTLIST_STRUCTURAL_COMPARATOR','control_definition':'all and only F3B3D final_best_tier=WEAK families from the same frozen 32-family priority panel; all formal family member loci included',
 'weak_family_n':9,'technical_unresolved_family_n':7,'technical_unresolved_treated_as_negative':False,
 'repeat_detector':'exactly same generic detector/parameters as F3B4B1: k=8,10,12,14,16; >=3 copies; start-gap 60-360; exact-core <=80; same strict/broad annotations',
 'threshold_retuning_after_survivor_results':False,'control_definition_frozen_before_control_results':True,'control_prespecified_before_b4b1_survivor_results':False,
 'f3b3d_summary_sha256':sha(d),'r4_family_map_sha256':sha(fm),'r4_instance_map_sha256':sha(im),'survivor_b4b1_family_summary_sha256':sha(sv),
 'implementation_sha256':sha(here/'scripts/weak_family_repeat_control.py'),
 'claim_boundary':'Post-discovery descriptive same-shortlist control. WEAK is not biochemical non-RT; unresolved remains censored.'}
(here/'frozen/F3B4B1C_IMPLEMENTATION_LOCK.json').write_text(json.dumps(lock,indent=2)+'\n')
print(json.dumps(lock,indent=2))
PY

echo '[2/2] scan all formal loci from the 9 WEAK families with unchanged B4B1 detector'
python "$HERE/scripts/weak_family_repeat_control.py" \
  "$D" "$FM" "$IM" "$S0/F3B1_orfs_aa" "$S0/F3B1_orfs_aa.index" "$S0/genomeDB" "$S0/genomeDB.index" \
  "$SV" "$HERE/outputs" "$HERE/frozen/F3B4B1C_IMPLEMENTATION_LOCK.json" 0 \
  2>&1 | tee "$HERE/logs/f3b4b1c_weak_control.log"
(cd "$HERE" && sha256sum frozen/F3B4B1C_IMPLEMENTATION_LOCK.json outputs/* > outputs/STAGE05F3B4B1C_SHA256.tsv)
echo PASS_STAGE05F3B4B1C_RUN_COMPLETE
