#!/usr/bin/env bash
set -euo pipefail

R5=${1:?R5 Stage05D4-R5 directory}
R4=${2:?R4 Stage05D4-R4 directory}
THREADS=${3:-32}

cd "$R5"

[ -s outputs/GATE_STAGE05D4R5_PREDICTION_LAYOUT.json ] || [ -s outputs/GATE_STAGE05D4R5_PREDICTION.json ] || {
  echo "[FAIL] no R5 prediction/layout gate found" >&2
  exit 2
}

TARGET="$R4/references/af_swissprot"
SPROT="$R4/references/uniprot_sprot.fasta.gz"
[ -f "${TARGET}.dbtype" ] || { echo "[FAIL] target DB missing: ${TARGET}.dbtype" >&2; exit 3; }
[ -s "$SPROT" ] || { echo "[FAIL] Swiss-Prot FASTA missing: $SPROT" >&2; exit 4; }
command -v foldseek >/dev/null 2>&1 || { echo "[FAIL] foldseek not in PATH" >&2; exit 5; }

mkdir -p work raw outputs logs tmp

# Preserve any prior downstream artifacts; do not delete scientific prediction files.
ts=$(date +%Y%m%d_%H%M%S)
for x in work/p1_cif_flat_regular work/P1_FULL3D_QUERY_DB work/P1_FULL3D_QUERY_DB.dbtype work/P1_FULL3D_QUERY_DB.index work/P1_FULL3D_QUERY_DB.lookup work/P1_FULL3D_QUERY_DB.source work/P1_FULL3D_QUERY_DB_h work/P1_FULL3D_QUERY_DB_h.dbtype work/P1_FULL3D_QUERY_DB_h.index work/P1_FULL3D_QUERY_DB_ca work/P1_FULL3D_QUERY_DB_ca.dbtype work/P1_FULL3D_QUERY_DB_ca.index work/P1_FULL3D_QUERY_DB_ss work/P1_FULL3D_QUERY_DB_ss.dbtype work/P1_FULL3D_QUERY_DB_ss.index work/P1_FULL3D_QUERY_DB.fasta raw/P1_FULL3D_AFSP.m8 tmp/full3d_search_v105; do
  if [ -e "$x" ] || [ -L "$x" ]; then
    mkdir -p "work/pre_v105_${ts}"
    bn=$(basename "$x")
    mv "$x" "work/pre_v105_${ts}/${bn}" 2>/dev/null || true
  fi
done

mkdir -p work/p1_cif_flat_regular

n_src=$(find boltz_out/predictions -type f -name '*_model_0.cif' | wc -l | tr -d ' ')
[ "$n_src" -eq 43 ] || { echo "[FAIL] expected 43 canonical regular CIFs, got $n_src" >&2; exit 6; }

# Build flat query directory with true regular files, not symlinks.
while IFS= read -r src; do
  q=$(basename "$(dirname "$src")")
  dst="work/p1_cif_flat_regular/${q}.cif"
  if ! ln "$src" "$dst" 2>/dev/null; then
    cp --reflink=auto "$src" "$dst"
  fi
done < <(find boltz_out/predictions -type f -name '*_model_0.cif' | sort)

n_flat=$(find work/p1_cif_flat_regular -maxdepth 1 -type f -name '*.cif' | wc -l | tr -d ' ')
n_link=$(find work/p1_cif_flat_regular -maxdepth 1 -type l | wc -l | tr -d ' ')
[ "$n_flat" -eq 43 ] || { echo "[FAIL] expected 43 flat regular CIFs, got $n_flat" >&2; exit 7; }
[ "$n_link" -eq 0 ] || { echo "[FAIL] flat query directory still contains $n_link symlinks" >&2; exit 8; }

echo "[PREFLIGHT] flat regular CIFs=$n_flat symlinks=$n_link"

# Explicitly create the Foldseek query DB. This avoids any ambiguity about
# directory traversal inside easy-search.
rm -rf tmp/querydb_v105
mkdir -p tmp/querydb_v105
foldseek createdb \
  work/p1_cif_flat_regular \
  work/P1_FULL3D_QUERY_DB \
  2>&1 | tee logs/foldseek_query_createdb_v105.log

[ -f work/P1_FULL3D_QUERY_DB.dbtype ] || { echo "[FAIL] query DB was not created" >&2; exit 9; }
[ -f work/P1_FULL3D_QUERY_DB_ca.dbtype ] || { echo "[FAIL] query C-alpha DB missing" >&2; exit 10; }
[ -f work/P1_FULL3D_QUERY_DB_ss.dbtype ] || { echo "[FAIL] query 3Di DB missing" >&2; exit 11; }

foldseek convert2fasta \
  work/P1_FULL3D_QUERY_DB \
  work/P1_FULL3D_QUERY_DB.fasta \
  2>&1 | tee logs/foldseek_query_convert2fasta_v105.log

n_db=$(grep -c '^>' work/P1_FULL3D_QUERY_DB.fasta || true)
[ "$n_db" -eq 43 ] || { echo "[FAIL] Foldseek query DB has $n_db entries, expected 43" >&2; exit 12; }
echo "[PASS] Foldseek query DB entries=$n_db"

rm -rf tmp/full3d_search_v105
mkdir -p tmp/full3d_search_v105
rm -f raw/P1_FULL3D_AFSP.m8

# Same frozen search settings as v1.0.0. Only the query input representation
# changed from an implicit directory to an explicit Foldseek structure DB.
foldseek easy-search \
  work/P1_FULL3D_QUERY_DB \
  "$TARGET" \
  raw/P1_FULL3D_AFSP.m8 \
  tmp/full3d_search_v105 \
  -e 0.01 \
  --max-seqs 50 \
  --threads "$THREADS" \
  --exact-tmscore 1 \
  --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,tmscore \
  2>&1 | tee logs/foldseek_full3d_v105.log

[ -e raw/P1_FULL3D_AFSP.m8 ] || { echo "[FAIL] Foldseek completed without creating raw/P1_FULL3D_AFSP.m8" >&2; exit 13; }

n_lines=$(wc -l < raw/P1_FULL3D_AFSP.m8 | tr -d ' ')
echo "[FOLDSEEK] m8_rows=$n_lines"

python scripts/summarize_r5.py \
  --freeze frozen \
  --boltz boltz_out \
  --m8 raw/P1_FULL3D_AFSP.m8 \
  --sprot "$SPROT" \
  --out outputs

python - <<'PY'
import json
p='outputs/GATE_STAGE05D4R5.json'
g=json.load(open(p))
assert g.get('gate_pass') is True, g
assert g.get('n_P1') == 43, g
assert g.get('n_boltz_predictions') == 43, g
print(json.dumps(g, indent=2))
print('[PASS] Stage05D4-R5 explicit-queryDB Foldseek + summarizer complete')
PY
