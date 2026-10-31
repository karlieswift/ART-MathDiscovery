#!/usr/bin/env bash
set -euo pipefail
R5=${1:?R5 package dir}
R4=${2:?R4 package dir}
THREADS=${3:-32}

cd "$R5"
mkdir -p raw outputs logs tmp work

QDIR=work/p1_cif_flat_regular
# Rebuild the flat query directory from the verified canonical regular-file layout.
rm -rf "$QDIR"
mkdir -p "$QDIR"
while IFS= read -r src; do
  base=$(basename "$src")
  dst="$QDIR/$base"
  ln "$src" "$dst" 2>/dev/null || cp -a "$src" "$dst"
done < <(find boltz_out/predictions -type f -name '*_model_0.cif' | sort)

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
TMP=tmp/full3d_search_v107
LOG=logs/foldseek_full3d_v107.log
rm -f "$OUTM8"
rm -rf "$TMP"
mkdir -p "$TMP"

# Foldseek 10.941cd33 does not implement a generic output token named 'tmscore'.
# Use the supported 'alntmscore' field. The existing summarizer reads column 11;
# therefore its legacy field name 'tmscore' means Foldseek alntmscore in v1.0.7.
foldseek easy-search \
  "$QDIR" \
  "$TARGET" \
  "$OUTM8" \
  "$TMP" \
  -e 0.01 \
  --max-seqs 50 \
  --threads "$THREADS" \
  --exact-tmscore 1 \
  --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore \
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
x['implementation_hotfix']='v1.0.7 downstream-only: replace unsupported Foldseek --format-output token tmscore with supported alntmscore; rebuild flat query directory from 43 verified canonical regular CIFs; no prediction rerun and no scientific threshold change.'
x['foldseek_tm_field']='alntmscore'
x['legacy_summary_field_note']='best_tmscore and best_rt_tmscore in v1.0.7 store Foldseek alntmscore because summarize_r5.py consumes the 11th m8 column.'
x['foldseek_query_representation']='43 regular mmCIF files rebuilt from boltz_out/predictions'
x['scientific_settings_changed']=False
x['prediction_workers_rerun']=False
json.dump(x,open(p,'w'),indent=2)
print(json.dumps(x,indent=2))
PY

python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/STAGE05D4R5_SUMMARY.json')
x=json.load(open(p))
x['foldseek_tm_field']='alntmscore'
x['legacy_summary_field_note']='best_tmscore and best_rt_tmscore store Foldseek alntmscore in v1.0.7.'
json.dump(x,open(p,'w'),indent=2)
PY

echo "[PASS] Stage05D4-R5 Foldseek search + summarizer complete with supported alntmscore field"
