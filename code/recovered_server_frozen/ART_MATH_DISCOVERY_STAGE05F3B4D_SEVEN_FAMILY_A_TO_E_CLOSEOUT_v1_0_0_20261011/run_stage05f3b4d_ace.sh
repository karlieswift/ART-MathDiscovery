#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
B4C="$P/ART_MATH_DISCOVERY_STAGE05F3B4C_GLOBAL_FORMAL_RT_ARCHITECTURE_CENSUS_AND_PARTNER_SCAN_v1_0_0_20261011/outputs"
LOCK="$S/frozen/F3B4D_IMPLEMENTATION_LOCK.json"
mkdir -p "$S"/{outputs,logs,work/known_reference,work/mmseqs_partner,work/mmseqs_known,work/mmseqs_global,work/mmseqs_art,raw}
for x in python mmseqs sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x"; exit 2; }; done
[[ -s "$LOCK" ]] || { echo '[FATAL] run prepare first'; exit 2; }

echo '[A-C/5] lineage dereplication + repeat-core conservation + base partner conservation'
python "$S/scripts/analyze_abc.py" \
  "$S/outputs/F3B4D_CORE_RECURRENT_STRICT_LOCI_16.tsv" \
  "$S/inputs/F3B4D_CORE_ARCHITECTURE_CONTEXT_16.fna" \
  "$B4C/F3B4C1_ALL_DETECTED_REPEAT_ARRAYS.tsv" \
  "$B4C/F3B4C2_RT_FAMILY_REPEAT_PARTNER_ARCHITECTURE_SUMMARY_4912.tsv" \
  "$S/inputs/F3B4D_CORE_PARTNER_16.faa" \
  "$S/outputs/F3B4D_TARGET_FAMILIES_7.tsv" \
  "$S/outputs" "$LOCK" \
  2>&1 | tee "$S/logs/abc_analysis.log"

echo '[C/5] descriptive partner catalog annotation'
PARTNER7="$S/inputs/F3B4D_PARTNER_REPRESENTATIVES_7.faa"
PFAM="$P/external_db/E3H_GLOBAL_CATALOG_20261008/pfam/Pfam-A.hmm"
PHROG="$P/external_db/E3H_GLOBAL_CATALOG_20261008/phrog/PHROG_all_proteins.faa"
VOG="$P/external_db/E3H_GLOBAL_CATALOG_20261008/vog/vogdb.proteins.all.faa"
PFOUT="$S/raw/F3B4D_C_PARTNER_VS_PFAM.tblout"
PHOUT="$S/raw/F3B4D_C_PARTNER_VS_PHROG.m8"
VGOUT="$S/raw/F3B4D_C_PARTNER_VS_VOG.m8"
: > "$PFOUT"; : > "$PHOUT"; : > "$VGOUT"
if command -v hmmscan >/dev/null 2>&1 && [[ -s "$PFAM" ]]; then
  set +e
  hmmscan --cpu "$THREADS" -E 1e-3 --domE 1e-3 --tblout "$PFOUT" "$PFAM" "$PARTNER7" > "$S/logs/partner_pfam.log" 2>&1
  rc=$?
  set -e
  [[ "$rc" -eq 0 ]] || { echo "[WARN] Pfam hmmscan technical failure rc=$rc"; : > "$PFOUT"; }
else echo '[WARN] Pfam resource/hmmscan unavailable; C remains sequence/topology complete, domain annotation partial'; fi
run_m8() {
  local target="$1" out="$2" tmp="$3" log="$4"
  if [[ -s "$target" ]]; then
    rm -rf "$tmp"; rm -f "$out"
    set +e
    mmseqs easy-search "$PARTNER7" "$target" "$out" "$tmp" --threads "$THREADS" -s 7.5 -e 1e-3 --max-seqs 20 --mask 0 --format-output 'query,target,fident,alnlen,qcov,tcov,evalue,bits' > "$log" 2>&1
    rc=$?
    set -e
    if [[ "$rc" -eq 0 && -s "$out" ]]; then { echo -e 'query\ttarget\tfident\talnlen\tqcov\ttcov\tevalue\tbits'; cat "$out"; } > "$out.tmp"; mv "$out.tmp" "$out"; else : > "$out"; echo "[WARN] catalog search failure target=$target rc=$rc"; fi
  fi
}
run_m8 "$PHROG" "$PHOUT" "$S/work/mmseqs_partner/phrog" "$S/logs/partner_phrog.log"
if [[ "${F3B4D_RUN_VOG:-0}" == 1 ]]; then run_m8 "$VOG" "$VGOUT" "$S/work/mmseqs_partner/vog" "$S/logs/partner_vog.log"; fi
python "$S/scripts/integrate_partner_catalogs.py" \
  "$S/outputs/F3B4D_ARCHITECTURE_REPRESENTATIVE_PER_FAMILY_7.tsv" \
  "$S/outputs/F3B4D_C_PARTNER_CONSERVATION_SUMMARY_7.tsv" \
  "$PFOUT" "$PHOUT" "$VGOUT" \
  "$S/outputs/F3B4D_C_PARTNER_CONSERVATION_AND_CATALOG_7.tsv" \
  "$S/outputs/STAGE05F3B4D_C_GATE.json" \
  2>&1 | tee "$S/logs/partner_catalog_summary.log"

echo '[E/5] frozen local known-system crosswalk; optional global/full-ART catalog if supplied'
python "$S/scripts/vendor_build_known_reference_panel.py" "$P" "$S/work/known_reference" > "$S/logs/e_reference_panel.log" 2>&1
REFN=$(python - "$S/work/known_reference/F3B3B_REFERENCE_PANEL_SUMMARY.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))['reference_unique_aa_n'])
PY
)
LOCAL="$S/raw/F3B4D_E_RT7_VS_LOCAL_KNOWN.m8"; : > "$LOCAL"
if [[ "$REFN" -gt 0 ]]; then
  rm -rf "$S/work/mmseqs_known/run"
  mmseqs easy-search "$S/inputs/F3B4D_RT_ARCHITECTURE_REPRESENTATIVES_7.faa" "$S/work/known_reference/F3B3B_KNOWN_SYSTEM_REFERENCE_PANEL.faa" "$LOCAL" "$S/work/mmseqs_known/run" --threads "$THREADS" -s 7.5 -e 1e-3 --max-seqs 100 --mask 0 --format-output 'query,target,fident,alnlen,qcov,tcov,evalue,bits' > "$S/logs/e_local_known.log" 2>&1
  if [[ -s "$LOCAL" ]]; then { echo -e 'query\ttarget\tfident\talnlen\tqcov\ttcov\tevalue\tbits'; cat "$LOCAL"; } > "$LOCAL.tmp"; mv "$LOCAL.tmp" "$LOCAL"; fi
fi
GLOBAL="$S/raw/F3B4D_E_RT7_VS_SUPPLIED_GLOBAL.m8"; ART="$S/raw/F3B4D_E_RT7_VS_SUPPLIED_ART.m8"; : > "$GLOBAL"; : > "$ART"
GP=0; AP=0
if [[ -n "${F3B4D_GLOBAL_REFERENCE_FASTAS:-}" ]]; then
  GP=1; GFA="$S/work/F3B4D_SUPPLIED_GLOBAL_REFERENCE.faa"; : > "$GFA"
  IFS=':' read -r -a GS <<< "$F3B4D_GLOBAL_REFERENCE_FASTAS"
  for f in "${GS[@]}"; do [[ -s "$f" ]] || { echo "[FATAL] missing supplied global FASTA $f"; exit 2; }; cat "$f" >> "$GFA"; done
  rm -rf "$S/work/mmseqs_global/run"
  mmseqs easy-search "$S/inputs/F3B4D_RT_ARCHITECTURE_REPRESENTATIVES_7.faa" "$GFA" "$GLOBAL" "$S/work/mmseqs_global/run" --threads "$THREADS" -s 7.5 -e 1e-3 --max-seqs 100 --mask 0 --format-output 'query,target,fident,alnlen,qcov,tcov,evalue,bits' > "$S/logs/e_global.log" 2>&1
  [[ ! -s "$GLOBAL" ]] || { { echo -e 'query\ttarget\tfident\talnlen\tqcov\ttcov\tevalue\tbits'; cat "$GLOBAL"; } > "$GLOBAL.tmp"; mv "$GLOBAL.tmp" "$GLOBAL"; }
fi
if [[ -n "${F3B4D_FULL_ART95_FASTA:-}" ]]; then
  AP=1; [[ -s "$F3B4D_FULL_ART95_FASTA" ]] || { echo '[FATAL] F3B4D_FULL_ART95_FASTA missing'; exit 2; }
  rm -rf "$S/work/mmseqs_art/run"
  mmseqs easy-search "$S/inputs/F3B4D_RT_ARCHITECTURE_REPRESENTATIVES_7.faa" "$F3B4D_FULL_ART95_FASTA" "$ART" "$S/work/mmseqs_art/run" --threads "$THREADS" -s 7.5 -e 1e-3 --max-seqs 100 --mask 0 --format-output 'query,target,fident,alnlen,qcov,tcov,evalue,bits' > "$S/logs/e_art95.log" 2>&1
  [[ ! -s "$ART" ]] || { { echo -e 'query\ttarget\tfident\talnlen\tqcov\ttcov\tevalue\tbits'; cat "$ART"; } > "$ART.tmp"; mv "$ART.tmp" "$ART"; }
fi
python "$S/scripts/summarize_e.py" \
  "$S/outputs/F3B4D_ARCHITECTURE_REPRESENTATIVE_PER_FAMILY_7.tsv" \
  "$S/work/known_reference/F3B3B_KNOWN_SYSTEM_REFERENCE_MAP.tsv" \
  "$LOCAL" "$GLOBAL" "$ART" "$GP" "$AP" \
  "$S/outputs/F3B4D_E_KNOWN_NOVELTY_CROSSWALK_7.tsv" "$S/outputs/STAGE05F3B4D_E_GATE.json" \
  2>&1 | tee "$S/logs/e_summary.log"

echo '[5/5] hashes'
sha256sum "$S/frozen/F3B4D_IMPLEMENTATION_LOCK.json" "$S/outputs/F3B4D_A_"* "$S/outputs/F3B4D_B_"* "$S/outputs/F3B4D_C_"* "$S/outputs/F3B4D_E_"* "$S/outputs/STAGE05F3B4D_A_GATE.json" "$S/outputs/STAGE05F3B4D_B_GATE.json" "$S/outputs/STAGE05F3B4D_C_GATE.json" "$S/outputs/STAGE05F3B4D_E_GATE.json" 2>/dev/null > "$S/outputs/STAGE05F3B4D_ACE_SHA256.tsv" || true
echo PASS_STAGE05F3B4D_ACE_RUN_COMPLETE
