#!/usr/bin/env bash
set -euo pipefail
# usage: ./run_stage05d4r5c.sh R5_ROOT [THREADS]
R5=${1:?R5_ROOT}; THREADS=${2:-32}
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
mkdir -p raw logs outputs tmp work/p1_cif_flat_regular
python - "$R5" <<'PY'
import json,sys,os
r=sys.argv[1];g=json.load(open(os.path.join(r,'outputs','GATE_STAGE05D4R5.json')))
assert g['gate_pass'] and g['n_P1']==43,g
f=json.load(open('frozen/STAGE05D4R5C_FREEZE.json'));assert f['gate_pass'],f
print('[PASS] R5 P1=43 and R5C freeze gate')
PY
rm -rf work/p1_cif_flat_regular; mkdir -p work/p1_cif_flat_regular
python - "$R5" <<'PY'
import os,glob,shutil,sys
r=sys.argv[1];files=sorted(glob.glob(os.path.join(r,'boltz_out','predictions','*','*_model_0.cif')))
if len(files)!=43: raise SystemExit(f'expected 43 CIFs got {len(files)}')
for p in files:
 q=os.path.basename(os.path.dirname(p)); dst=os.path.join('work','p1_cif_flat_regular',q+'_model_0.cif')
 try: os.link(p,dst)
 except OSError: shutil.copy2(p,dst)
print('[PASS] staged 43 P1 CIFs')
PY
# P1 vs clean annotation-confirmed RT structures.
rm -rf tmp/p1_clean_rt
foldseek easy-search work/p1_cif_flat_regular references/clean_rt_afsp raw/P1_CLEAN_RT.m8 tmp/p1_clean_rt -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore 2>&1 | tee logs/foldseek_p1_clean_rt.log
# DB-vs-DB calibration searches. Search returns DB, convertalis emits requested columns.
rm -rf tmp/pos_search tmp/neg_search raw/POS_VS_CLEAN.db* raw/NEG_VS_CLEAN.db*
foldseek search references/clean_rt_afsp references/clean_rt_afsp raw/POS_VS_CLEAN.db tmp/pos_search -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 -a > logs/foldseek_pos_search.log 2>&1
foldseek convertalis references/clean_rt_afsp references/clean_rt_afsp raw/POS_VS_CLEAN.db raw/POS_VS_CLEAN.m8 --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore --exact-tmscore 1 > logs/foldseek_pos_convertalis.log 2>&1
foldseek search references/matched_nonrt_afsp references/clean_rt_afsp raw/NEG_VS_CLEAN.db tmp/neg_search -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 -a > logs/foldseek_neg_search.log 2>&1
foldseek convertalis references/matched_nonrt_afsp references/clean_rt_afsp raw/NEG_VS_CLEAN.db raw/NEG_VS_CLEAN.m8 --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore --exact-tmscore 1 > logs/foldseek_neg_convertalis.log 2>&1
printf '[ROWS] P1=%s POS=%s NEG=%s\n' "$(wc -l < raw/P1_CLEAN_RT.m8)" "$(wc -l < raw/POS_VS_CLEAN.m8)" "$(wc -l < raw/NEG_VS_CLEAN.m8)"
python scripts/summarize_r5c.py --p1m8 raw/P1_CLEAN_RT.m8 --posm8 raw/POS_VS_CLEAN.m8 --negm8 raw/NEG_VS_CLEAN.m8 --p1freeze "$R5/frozen/STAGE05D4R5_P1_FREEZE.tsv" --pos frozen/STAGE05D4R5C_CLEAN_RT_REFERENCES.tsv --neg frozen/STAGE05D4R5C_MATCHED_NONRT_CONTROLS.tsv --outdir outputs
sha256sum raw/P1_CLEAN_RT.m8 raw/POS_VS_CLEAN.m8 raw/NEG_VS_CLEAN.m8 outputs/STAGE05D4R5C_*.tsv outputs/STAGE05D4R5C_SUMMARY.json outputs/GATE_STAGE05D4R5C.json > outputs/STAGE05D4R5C_SHA256.tsv
cat outputs/STAGE05D4R5C_SUMMARY.json
