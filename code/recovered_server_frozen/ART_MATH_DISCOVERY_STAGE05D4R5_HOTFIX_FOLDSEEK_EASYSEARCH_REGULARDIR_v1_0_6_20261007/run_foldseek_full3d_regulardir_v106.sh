#!/usr/bin/env bash
set -euo pipefail
R5=${1:?R5 package dir}
R4=${2:?R4 package dir}
THREADS=${3:-32}

cd "$R5"
mkdir -p raw outputs logs tmp work

QDIR=work/p1_cif_flat_regular
[ -d "$QDIR" ] || { echo "[FAIL] missing $R5/$QDIR; run v1.0.5 preflight/createdb once first"; exit 2; }

n_cif=$(find "$QDIR" -maxdepth 1 -type f -name '*.cif' | wc -l)
n_link=$(find "$QDIR" -maxdepth 1 -type l | wc -l)
echo "[PREFLIGHT] cwd=$PWD"
echo "[PREFLIGHT] regular_CIF=$n_cif symlinks=$n_link"
[ "$n_cif" -eq 43 ] || { echo "[FAIL] expected 43 regular CIFs, got $n_cif"; exit 3; }
[ "$n_link" -eq 0 ] || { echo "[FAIL] expected 0 symlinks, got $n_link"; exit 4; }

TARGET="$R4/references/af_swissprot"
[ -e "$TARGET" ] || [ -e "${TARGET}.dbtype" ] || { echo "[FAIL] missing target Foldseek DB prefix: $TARGET"; exit 5; }
SPROT="$R4/references/uniprot_sprot.fasta.gz"
[ -f "$SPROT" ] || { echo "[FAIL] missing $SPROT"; exit 6; }

OUTM8=raw/P1_FULL3D_AFSP.m8
TMP=tmp/full3d_search_v106
LOG=logs/foldseek_full3d_v106.log
rm -f "$OUTM8"
rm -rf "$TMP"
mkdir -p "$TMP"

# IMPORTANT: easy-search accepts raw structure files/folders as query input and a prebuilt DB as target.
# Do NOT pass the prebuilt queryDB prefix to easy-search; DB-vs-DB uses foldseek search instead.
foldseek easy-search \
  "$QDIR" \
  "$TARGET" \
  "$OUTM8" \
  "$TMP" \
  -e 0.01 \
  --max-seqs 50 \
  --threads "$THREADS" \
  --exact-tmscore 1 \
  --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,tmscore \
  2>&1 | tee "$LOG"

[ -f "$OUTM8" ] || { echo "[FAIL] Foldseek returned without creating $OUTM8"; exit 7; }
rows=$(wc -l < "$OUTM8")
echo "[FOLDSEEK] m8_rows=$rows"

python scripts/summarize_r5.py \
  --freeze frozen \
  --boltz boltz_out \
  --m8 "$OUTM8" \
  --sprot "$SPROT" \
  --out outputs

python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/GATE_STAGE05D4R5.json')
x=json.load(open(p))
x['implementation_hotfix']='v1.0.6 downstream-only: run Foldseek easy-search on the verified directory of 43 regular CIF query files; v1.0.5 queryDB creation remains audit-only because easy-search expects raw query files/folders, while DB-vs-DB uses foldseek search.'
x['foldseek_query_representation']='43 regular mmCIF files in work/p1_cif_flat_regular'
x['scientific_settings_changed']=False
x['prediction_workers_rerun']=False
json.dump(x,open(p,'w'),indent=2)
print(json.dumps(x,indent=2))
PY

echo "[PASS] Stage05D4-R5 Foldseek regular-directory search + summarizer complete"
