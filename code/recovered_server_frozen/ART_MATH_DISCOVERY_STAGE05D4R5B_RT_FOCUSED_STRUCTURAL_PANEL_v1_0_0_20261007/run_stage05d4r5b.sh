#!/usr/bin/env bash
set -euo pipefail
# usage: ./run_stage05d4r5b.sh R5_ROOT [THREADS]
R5=${1:?R5_ROOT}
THREADS=${2:-32}
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
mkdir -p outputs raw logs tmp/search work/p1_cif_flat_regular
python - "$R5" <<'PY'
import json,sys,os
r=sys.argv[1]
g=json.load(open(os.path.join(r,'outputs','GATE_STAGE05D4R5.json')))
assert g['gate_pass'],g
assert g['n_P1']==43,g
print('[PASS] R5 final gate and P1=43')
f=json.load(open('frozen/STAGE05D4R5B_REFERENCE_FREEZE.json'))
assert f['gate_pass'] and f['gate']=='PASS_RT_FOCUSED_REFERENCE_PANEL_FROZEN',f
assert f['reference_panel_selected_from_candidate_results'] is False,f
print('[PASS] R5B external reference panel frozen before search')
PY
rm -rf work/p1_cif_flat_regular
mkdir -p work/p1_cif_flat_regular
python - "$R5" <<'PY'
import os,glob,shutil,sys
r=sys.argv[1]
files=sorted(glob.glob(os.path.join(r,'boltz_out','predictions','*','*_model_0.cif')))
if len(files)!=43: raise SystemExit(f'expected 43 R5 CIFs, got {len(files)}')
for p in files:
    q=os.path.basename(os.path.dirname(p)); dst=os.path.join('work','p1_cif_flat_regular',q+'_model_0.cif')
    try: os.link(p,dst)
    except OSError: shutil.copy2(p,dst)
print('[PASS] staged 43 regular full-coordinate query CIFs')
PY
rm -rf tmp/search
foldseek easy-search work/p1_cif_flat_regular references/rt_panel_afsp raw/P1_RT_PANEL.m8 tmp/search \
  -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 \
  --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore \
  2>&1 | tee logs/foldseek_r5b.log
python scripts/summarize_r5b.py \
  --m8 raw/P1_RT_PANEL.m8 \
  --freeze "$R5/frozen/STAGE05D4R5_P1_FREEZE.tsv" \
  --panel frozen/STAGE05D4R5B_RT_REFERENCE_PANEL.tsv \
  --outdir outputs
sha256sum raw/P1_RT_PANEL.m8 outputs/STAGE05D4R5B_P1_RT_STRUCTURAL_TRIAGE.tsv outputs/STAGE05D4R5B_ALL_RT_PANEL_HITS.tsv outputs/STAGE05D4R5B_SUMMARY.json outputs/GATE_STAGE05D4R5B.json > outputs/STAGE05D4R5B_SHA256.tsv
cat outputs/STAGE05D4R5B_SUMMARY.json
