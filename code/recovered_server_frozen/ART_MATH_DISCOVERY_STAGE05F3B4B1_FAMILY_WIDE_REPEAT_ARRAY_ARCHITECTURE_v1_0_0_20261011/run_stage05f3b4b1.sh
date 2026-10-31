#!/usr/bin/env bash
set -euo pipefail
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/outputs" "$HERE/logs" "$HERE/frozen"
B0="$P/ART_MATH_DISCOVERY_STAGE05F3B4B0_FAMILY_WIDE_CONTEXT_COMPLETENESS_RESCUE_v1_0_0_20261011/outputs"
ALL="$B0/F3B4B0_ALL_SURVIVOR_FAMILY_MEMBER_LOCUS_COMPLETENESS.tsv"
FS="$B0/F3B4B0_SURVIVOR_FAMILY_CONTEXT_COMPLETENESS_SUMMARY.tsv"
GATE="$B0/STAGE05F3B4B0_GATE.json"
S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010/work"
for f in "$ALL" "$FS" "$GATE" "$S0/genomeDB" "$S0/genomeDB.index"; do [[ -e "$f" ]] || { echo "[FATAL] missing $f"; exit 2; }; done

"$HERE/tests/run_selftest.sh"

echo '[1/3] freeze prospective repeat-array scan before any repeat result'
python - "$HERE" "$ALL" "$FS" "$GATE" <<'PY'
import json,hashlib,sys
from pathlib import Path
here,allp,fsp,gp=map(Path,sys.argv[1:])
def sha(p):
 h=hashlib.sha256(); h.update(p.read_bytes()); return h.hexdigest()
g=json.load(open(gp)); assert g['gate_pass'] is True and g['selected_family_candidate_instance_n']==1605
assert sha(allp)=='acc14bf1b0f8f7ca44769e93051ff4743bfa63db50618ca9fed14360e1ea7692'
assert sha(fsp)=='5edf27b1fffe5c5a1b11e82e7e3af5aedc688179535eec088b5b8d270fa2c4f6'
lock={
 'stage':'Stage05F3B4B1-implementation-freeze','gate':'PASS_F3B4B1_REPEAT_ARRAY_SCAN_FROZEN_PRE_REPEAT_RESULTS','gate_pass':True,
 'input_family_n':16,'input_locus_n':1605,
 'input_rule':'all and only 1605 formal candidate instances in the 16 frozen F3B3D STRONG/SUGGESTIVE survivor families; no repeat-based locus preselection',
 'candidate_orientation':'all sequences transformed to candidate 5-prime-to-3-prime orientation before UPSTREAM/DOWNSTREAM labeling',
 'scan_flank_each_side_bp':10000,'seed_k_values':[8,10,12,14,16],'min_repeat_copies':3,
 'allowed_start_to_start_gap_bp':[60,360],'max_exact_core_motif_bp':80,
 'generic_array_role':'primary descriptive discovery screen; no Anthropic geometry required for detection',
 'strict_published_art_geometry_annotation':'UPSTREAM; 3-21 copies; every inferred spacer 120-220 bp',
 'broad_art_like_geometry_annotation':'UPSTREAM; 3-30 copies; median spacer 100-240 bp and >=75% spacers 80-280 bp',
 'published_art_geometry_used_for_candidate_filtering':False,
 'family_best_array_rule':'copy_n descending, motif_len descending, spacer_cv ascending, array_span descending, upstream before downstream, motif/candidate lexical tie-break',
 'absence_call_rule':'presence may be reported from partial contexts; no family-level absence claim; per-locus no-array observation is only labeled when >=5 kb upstream context is present',
 'partner_analysis_performed':False,'known_art_membership_used':False,'novelty_used':False,
 'b4b0_all_instances_sha256':sha(allp),'b4b0_family_summary_sha256':sha(fsp),'b4b0_gate_sha256':sha(gp),
 'implementation_sha256':sha(here/'scripts/repeat_array_scan.py'),
 'claim_boundary':'Frozen repeat-array architecture discovery only. ART-like geometry is an annotation, not a definition of novelty or function.'}
(here/'frozen/F3B4B1_IMPLEMENTATION_LOCK.json').write_text(json.dumps(lock,indent=2)+'\n')
print(json.dumps(lock,indent=2))
PY

echo '[2/3] scan repeat arrays across all 1605 exact loci'
python "$HERE/scripts/repeat_array_scan.py" \
  "$ALL" "$S0/genomeDB" "$S0/genomeDB.index" \
  "$HERE/outputs" "$HERE/frozen/F3B4B1_IMPLEMENTATION_LOCK.json" "$FS" \
  2>&1 | tee "$HERE/logs/f3b4b1_repeat_array_scan.log"

echo '[3/3] hashes'
(cd "$HERE" && sha256sum frozen/F3B4B1_IMPLEMENTATION_LOCK.json outputs/* > outputs/STAGE05F3B4B1_SHA256.tsv)
echo PASS_STAGE05F3B4B1_RUN_COMPLETE
