#!/usr/bin/env bash
set -euo pipefail
# usage: ./freeze_reference_panel.sh PROJECT_ROOT R4_ROOT RVT_ALL_HMM [THREADS]
P=${1:?PROJECT_ROOT}
R4=${2:?R4_ROOT}
HMM=${3:?RVT_ALL_HMM}
THREADS=${4:-32}
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
mkdir -p frozen raw logs work references tmp
command -v hmmsearch >/dev/null || { echo '[ERROR] hmmsearch not found'; exit 2; }
command -v foldseek >/dev/null || { echo '[ERROR] foldseek not found'; exit 2; }
[ -s "$HMM" ] || { echo "[ERROR] missing RVT-All HMM: $HMM"; exit 2; }
SP="$R4/references/uniprot_sprot.fasta.gz"
AF="$R4/references/af_swissprot"
[ -s "$SP" ] || { echo "[ERROR] missing frozen Swiss-Prot FASTA: $SP"; exit 2; }
[ -e "$AF" ] || { echo "[ERROR] missing frozen AF/Swiss-Prot DB: $AF"; exit 2; }
[ -s "$AF.lookup" ] || { echo "[ERROR] missing AF/Swiss-Prot lookup: $AF.lookup"; exit 2; }
EH='74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e'
ES='a9c3496aa727e61d25d9e1291c9ba28d2c2ad19c0e26569dc083d53ec4169536'
AH=$(sha256sum "$HMM" | awk '{print $1}')
AS=$(sha256sum "$SP" | awk '{print $1}')
[ "$AH" = "$EH" ] || { echo "[ERROR] RVT-All hash mismatch $AH"; exit 2; }
[ "$AS" = "$ES" ] || { echo "[ERROR] Swiss-Prot hash mismatch $AS"; exit 2; }
cp -f "$HMM" references/RVT-All.hmm
ln -sfn "$SP" references/uniprot_sprot.fasta.gz
ln -sfn "$AF" references/af_swissprot
ln -sfn "$AF.lookup" references/af_swissprot.lookup
printf '%s  %s\n' "$AH" references/RVT-All.hmm > frozen/SOURCE_RVT_ALL.sha256
printf '%s  %s\n' "$AS" references/uniprot_sprot.fasta.gz > frozen/SOURCE_SWISSPROT.sha256
sha256sum "$AF.lookup" > frozen/SOURCE_AFSP_LOOKUP.sha256

gzip -cd "$SP" > work/uniprot_sprot.fasta
# External-reference screen only; no P1 candidate sequence is used here.
hmmsearch --cpu "$THREADS" --noali -E 1e-5 --domE 1e-5 \
  --tblout raw/RVTALL_SWISSPROT.tblout \
  --domtblout raw/RVTALL_SWISSPROT.domtblout \
  references/RVT-All.hmm work/uniprot_sprot.fasta \
  > logs/hmmsearch_rvtall_swissprot.log 2>&1
python scripts/build_rt_panel.py \
  --sprot "$SP" \
  --hmmtbl raw/RVTALL_SWISSPROT.tblout \
  --lookup "$AF.lookup" \
  --top-per-profile 10 \
  --outdir frozen
python - <<'PY'
import json
s=json.load(open('frozen/STAGE05D4R5B_REFERENCE_FREEZE.json'))
assert s['gate_pass'],s
assert s['reference_panel_selected_from_candidate_results'] is False,s
print('[PASS] reference panel freeze:',s)
PY
rm -rf references/rt_panel_afsp*
foldseek createsubdb frozen/STAGE05D4R5B_RT_REFERENCE_KEYS.txt "$AF" references/rt_panel_afsp --subdb-mode 1 \
  > logs/foldseek_createsubdb.log 2>&1
[ -s references/rt_panel_afsp.dbtype ] || { echo '[ERROR] subset DB missing'; exit 2; }
[ -s references/rt_panel_afsp_ss.dbtype ] || { echo '[ERROR] subset 3Di DB missing'; exit 2; }
[ -s references/rt_panel_afsp_ca.dbtype ] || { echo '[ERROR] subset C-alpha DB missing'; exit 2; }
foldseek convert2fasta references/rt_panel_afsp work/rt_panel.faa > logs/foldseek_panel_convert2fasta.log 2>&1
N=$(grep -c '^>' work/rt_panel.faa || true)
E=$(python - <<'PY'
import json
print(json.load(open('frozen/STAGE05D4R5B_REFERENCE_FREEZE.json'))['afsp_structural_reference_entries'])
PY
)
[ "$N" = "$E" ] || { echo "[ERROR] panel DB count mismatch expected=$E got=$N"; exit 2; }
sha256sum frozen/STAGE05D4R5B_RT_REFERENCE_PANEL.tsv frozen/STAGE05D4R5B_RT_REFERENCE_KEYS.txt frozen/STAGE05D4R5B_PROFILE_COVERAGE.tsv frozen/STAGE05D4R5B_REFERENCE_FREEZE.json > frozen/STAGE05D4R5B_REFERENCE_FREEZE.sha256
printf '[PASS] RT-focused reference DB frozen entries=%s\n' "$N"
