#!/usr/bin/env bash
set -euo pipefail
# usage: ./freeze_stage05d4r5c.sh R4_ROOT R5B_ROOT [THREADS]
R4=${1:?R4_ROOT}; R5B=${2:?R5B_ROOT}; THREADS=${3:-32}
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
mkdir -p frozen references logs work tmp
python - "$R5B" <<'PY'
import json,sys,os
r=sys.argv[1]
g=json.load(open(os.path.join(r,'outputs','GATE_STAGE05D4R5B.json')))
assert g['gate_pass'] and g['gate']=='PASS_P1_RT_FOCUSED_STRUCTURAL_PANEL_SEARCH_COMPLETE',g
f=json.load(open(os.path.join(r,'frozen','STAGE05D4R5B_REFERENCE_FREEZE.json')))
assert f['gate_pass'] and f['reference_panel_selected_from_candidate_results'] is False,f
print('[PASS] source R5B frozen panel/gate')
PY
SP="$R4/references/uniprot_sprot.fasta.gz"; AF="$R4/references/af_swissprot"; PANEL="$R5B/frozen/STAGE05D4R5B_RT_REFERENCE_PANEL.tsv"
[ -s "$SP" ] || { echo '[ERROR] missing Swiss-Prot'; exit 2; }
[ -e "$AF" ] || { echo '[ERROR] missing AF/Swiss-Prot DB'; exit 2; }
[ -s "$AF.lookup" ] || { echo '[ERROR] missing AF lookup'; exit 2; }
[ -s "$PANEL" ] || { echo '[ERROR] missing R5B panel'; exit 2; }
python scripts/build_calibration_panel.py --panel "$PANEL" --sprot "$SP" --lookup "$AF.lookup" --outdir frozen --per-rt 3 --seed ART_STAGE05D4R5C_20261007
python - <<'PY'
import json
s=json.load(open('frozen/STAGE05D4R5C_FREEZE.json'));assert s['gate_pass'],s
print('[PASS] calibration universe freeze:',s)
PY
rm -rf references/clean_rt_afsp* references/matched_nonrt_afsp*
foldseek createsubdb frozen/STAGE05D4R5C_CLEAN_RT_KEYS.txt "$AF" references/clean_rt_afsp --subdb-mode 1 > logs/createsubdb_clean_rt.log 2>&1
foldseek createsubdb frozen/STAGE05D4R5C_NONRT_KEYS.txt "$AF" references/matched_nonrt_afsp --subdb-mode 1 > logs/createsubdb_nonrt.log 2>&1
for D in references/clean_rt_afsp references/matched_nonrt_afsp; do
  [ -s "$D.dbtype" ] && [ -s "${D}_ss.dbtype" ] && [ -s "${D}_ca.dbtype" ] || { echo "[ERROR] incomplete DB $D"; exit 2; }
done
foldseek convert2fasta references/clean_rt_afsp work/clean_rt.faa > logs/convert_clean_rt.log 2>&1
foldseek convert2fasta references/matched_nonrt_afsp work/nonrt.faa > logs/convert_nonrt.log 2>&1
NRT=$(grep -c '^>' work/clean_rt.faa || true); NNEG=$(grep -c '^>' work/nonrt.faa || true)
EXP_RT=$(python -c "import json;print(json.load(open('frozen/STAGE05D4R5C_FREEZE.json'))['n_clean_rt_structures'])")
EXP_NEG=$(python -c "import json;print(json.load(open('frozen/STAGE05D4R5C_FREEZE.json'))['n_matched_nonrt_controls'])")
[ "$NRT" = "$EXP_RT" ] || { echo "[ERROR] clean RT DB count $NRT != $EXP_RT"; exit 2; }
[ "$NNEG" = "$EXP_NEG" ] || { echo "[ERROR] nonRT DB count $NNEG != $EXP_NEG"; exit 2; }
sha256sum frozen/STAGE05D4R5C_* > frozen/STAGE05D4R5C_FREEZE_SHA256.tsv
printf '[PASS] Stage05D4-R5C frozen cleanRT=%s matchedNonRT=%s\n' "$NRT" "$NNEG"
