#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
B4C="$P/ART_MATH_DISCOVERY_STAGE05F3B4C_GLOBAL_FORMAL_RT_ARCHITECTURE_CENSUS_AND_PARTNER_SCAN_v1_0_0_20261011/outputs"
R4="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_v1_0_0_20261011/outputs"
S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010/work"
ER1="$P/ART_MATH_DISCOVERY_STAGE05F3B3ER1_FOLDSEEK_ONLY_TECHNICAL_REPAIR_v1_0_0_20261011/outputs/F3B3ER1_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv"
LOCK="$S/frozen/F3B4D_IMPLEMENTATION_LOCK.json"
mkdir -p "$S"/{frozen,outputs,logs,inputs/structure_yaml,work,raw/foldseek_d,tmp/foldseek_d,state,structure_boltz_out}
for x in python sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x"; exit 2; }; done
for f in \
  "$B4C/STAGE05F3B4C2_GATE.json" \
  "$B4C/F3B4C2_RT_FAMILY_REPEAT_PARTNER_ARCHITECTURE_SUMMARY_4912.tsv" \
  "$B4C/F3B4C2_RT_REPEAT_PARTNER_LOCUS_ARCHITECTURE_10598.tsv" \
  "$B4C/F3B4C1_ALL_DETECTED_REPEAT_ARRAYS.tsv" \
  "$B4C/F3B4C1_GLOBAL_FORMAL_RT_LOCUS_REPEAT_CENSUS_10598.tsv" \
  "$B4C/F3B4C2_NEAREST_DOWNSTREAM_ORF_PER_RT_LOCUS.tsv" \
  "$B4C/F3B4C2_NEAREST_DOWNSTREAM_ORF_PROTEINS.faa" \
  "$R4/F3B1R3_INSTANCE_MAP_WITH_CANONICAL_SOURCE.tsv" \
  "$S0/F3B1_orfs_aa" "$S0/F3B1_orfs_aa.index" "$S0/genomeDB" "$S0/genomeDB.index" \
  "$ER1" "$S/scripts/vendor_context_helper.py" "$S/scripts/prepare_targets.py"
do [[ -s "$f" ]] || { echo "[FATAL] missing $f"; exit 2; }; done

if [[ ! -s "$LOCK" ]]; then
  echo '[1/2] freeze seven-family A-E target/rules before A-E closeout results'
  python "$S/scripts/prepare_targets.py" \
    "$B4C" \
    "$R4/F3B1R3_INSTANCE_MAP_WITH_CANONICAL_SOURCE.tsv" \
    "$S0/F3B1_orfs_aa" "$S0/F3B1_orfs_aa.index" \
    "$S0/genomeDB" "$S0/genomeDB.index" \
    "$ER1" "$S/scripts/vendor_context_helper.py" \
    "$S/outputs" "$S/inputs" "$S/frozen" "$S/inputs/structure_yaml" "$LOCK" "$P" \
    2>&1 | tee "$S/logs/prepare_and_freeze.log"
  sha256sum "$LOCK" "$S/scripts"/*.py "$S/run_stage05f3b4d_prepare.sh" > "$S/frozen/F3B4D_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$S/outputs/F3B4D_TARGET_FAMILIES_7.tsv" "$S/outputs/F3B4D_CORE_RECURRENT_STRICT_LOCI_16.tsv" "$S/outputs/F3B4D_D_STRUCTURE_EXTENSION_TARGETS_6.tsv"
else
  echo '[1/2] [SKIP] implementation lock exists'
fi

echo '[2/2] compact manifests'
cat "$S/outputs/F3B4D_TARGET_FAMILIES_7.tsv"
echo
cat "$S/outputs/F3B4D_D_STRUCTURE_EXTENSION_TARGETS_6.tsv"
echo PASS_STAGE05F3B4D_PREPARE_RUN_COMPLETE
