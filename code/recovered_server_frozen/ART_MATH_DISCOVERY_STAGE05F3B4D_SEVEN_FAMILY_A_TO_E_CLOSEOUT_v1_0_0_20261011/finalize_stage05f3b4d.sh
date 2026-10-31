#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
for f in \
 "$S/outputs/F3B4D_TARGET_FAMILIES_7.tsv" \
 "$S/outputs/F3B4D_A_FAMILY_CONTEXT_DEREP_SUMMARY_7.tsv" \
 "$S/outputs/F3B4D_B_REPEAT_CORE_CONSERVATION_SUMMARY_7.tsv" \
 "$S/outputs/F3B4D_C_PARTNER_CONSERVATION_AND_CATALOG_7.tsv" \
 "$S/outputs/F3B4D_D_STRUCTURE_R5C_SUMMARY_6.tsv" \
 "$S/outputs/F3B4D_E_KNOWN_NOVELTY_CROSSWALK_7.tsv" \
 "$S/outputs/STAGE05F3B4D_A_GATE.json" "$S/outputs/STAGE05F3B4D_B_GATE.json" "$S/outputs/STAGE05F3B4D_C_GATE.json" "$S/outputs/STAGE05F3B4D_D_GATE.json" "$S/outputs/STAGE05F3B4D_E_GATE.json"
do [[ -s "$f" ]] || { echo "[FATAL] missing $f"; exit 2; }; done
python "$S/scripts/finalize_closeout.py" \
 "$S/outputs/F3B4D_TARGET_FAMILIES_7.tsv" \
 "$S/outputs/F3B4D_A_FAMILY_CONTEXT_DEREP_SUMMARY_7.tsv" \
 "$S/outputs/F3B4D_B_REPEAT_CORE_CONSERVATION_SUMMARY_7.tsv" \
 "$S/outputs/F3B4D_C_PARTNER_CONSERVATION_AND_CATALOG_7.tsv" \
 "$S/outputs/F3B4D_D_STRUCTURE_R5C_SUMMARY_6.tsv" \
 "$S/outputs/F3B4D_E_KNOWN_NOVELTY_CROSSWALK_7.tsv" \
 "$S/outputs/STAGE05F3B4D_A_GATE.json" "$S/outputs/STAGE05F3B4D_B_GATE.json" "$S/outputs/STAGE05F3B4D_C_GATE.json" "$S/outputs/STAGE05F3B4D_D_GATE.json" "$S/outputs/STAGE05F3B4D_E_GATE.json" \
 "$S/outputs" | tee "$S/logs/final_closeout.log"
sha256sum "$S/outputs/F3B4D_FINAL_SEVEN_FAMILY_COMPUTATIONAL_EVIDENCE_BOARD.tsv" "$S/outputs/STAGE05F3B4D_GATE.json" > "$S/outputs/STAGE05F3B4D_FINAL_SHA256.tsv"
echo PASS_STAGE05F3B4D_FINALIZE_COMPLETE
