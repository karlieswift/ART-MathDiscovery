#!/usr/bin/env bash
set -euo pipefail
# usage: ./run_stage05f2d.sh F1_INPUT F2A_INPUT F2CB_INPUT_OR_NONE [THREADS]
F1=${1:?F1 input required}
F2A=${2:?F2A input required}
F2CB=${3:-NONE}
THREADS=${4:-32}
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
rm -rf outputs raw logs work
mkdir -p outputs raw logs work/mmseqs_tmp
python scripts/stage05f2d.py prepare --f1 "$F1" --f2a "$F2A" --f2cb "$F2CB" --outdir outputs --work work | tee logs/prepare.log
command -v mmseqs >/dev/null 2>&1 || { echo '[FAIL] mmseqs not found in PATH' >&2; exit 2; }
python - <<'PY'
import json
p=json.load(open('work/RESOLVED_INPUTS.json'))
print(p['query_fasta'])
print(p['target_fasta'])
PY
QUERY=$(python - <<'PY'
import json; print(json.load(open('work/RESOLVED_INPUTS.json'))['query_fasta'])
PY
)
TARGET=$(python - <<'PY'
import json; print(json.load(open('work/RESOLVED_INPUTS.json'))['target_fasta'])
PY
)
mmseqs easy-search \
  "$QUERY" \
  "$TARGET" \
  raw/STAGE05F2D_QUERY30_VS_F1_UNIQUE.m8 \
  work/mmseqs_tmp \
  --threads "$THREADS" \
  -s 7.5 \
  -e 0.01 \
  --max-seqs 10000 \
  --format-output query,target,evalue,bits,fident,alnlen,qcov,tcov,qlen,tlen \
  2>&1 | tee logs/mmseqs_family_expansion.log
python scripts/stage05f2d.py summarize --m8 raw/STAGE05F2D_QUERY30_VS_F1_UNIQUE.m8 --outdir outputs --work work | tee logs/summarize.log
cat outputs/STAGE05F2D_GATE.json
