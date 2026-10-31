#!/usr/bin/env bash
set -euo pipefail
R5C=${1:?R5C_ROOT}; THREADS=${2:-32}
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
python - "$R5C" <<'PY'
import json,os,sys
r=sys.argv[1];g=json.load(open(os.path.join(r,'outputs','GATE_STAGE05D4R5C.json')))
assert g['gate_pass'] and g['calibration_informative'] is True,g
assert g['known_rt_suggestive_or_strong_sensitivity']>=0.75 and g['matched_nonrt_suggestive_or_strong_rate']<=0.05,g
assert os.path.exists(os.path.join(r,'references','clean_rt_afsp.dbtype')) or os.path.exists(os.path.join(r,'references','clean_rt_afsp'))
print('[PASS] informative R5C calibration and clean-RT reference DB available')
PY
[ -f outputs/GATE_STAGE05D4R5D_PREDICTION.json ] || { echo '[FAIL] P2 prediction gate missing'; exit 2; }
python - <<'PY'
import json
x=json.load(open('outputs/GATE_STAGE05D4R5D_PREDICTION.json'));assert x['gate_pass'] and x['n_P2']==40 and x['n_cif']==40 and x['n_confidence_json']==40,x
PY
rm -rf work/p2_cif_flat_regular tmp/p2_clean_rt
mkdir -p work/p2_cif_flat_regular raw logs
python - <<'PY'
import glob,os,shutil
fs=sorted(glob.glob('boltz_p2_out/predictions/*/*_model_0.cif'))
assert len(fs)==40,len(fs)
for p in fs:
 q=os.path.basename(os.path.dirname(p));dst=os.path.join('work','p2_cif_flat_regular',q+'_model_0.cif')
 try: os.link(p,dst)
 except OSError: shutil.copy2(p,dst)
print('[PASS] staged 40 P2 regular CIFs')
PY
foldseek easy-search work/p2_cif_flat_regular "$R5C/references/clean_rt_afsp" raw/P2_CLEAN_RT.m8 tmp/p2_clean_rt -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore 2>&1 | tee logs/foldseek_p2_clean_rt.log
printf '[ROWS] P2_CLEAN_RT=%s\n' "$(wc -l < raw/P2_CLEAN_RT.m8)"
python scripts/summarize_p2.py --freeze frozen/STAGE05D4R5D_P2_FREEZE.tsv --boltz boltz_p2_out --m8 raw/P2_CLEAN_RT.m8 --r5csummary "$R5C/outputs/STAGE05D4R5C_SUMMARY.json" --outdir outputs
sha256sum raw/P2_CLEAN_RT.m8 outputs/STAGE05D4R5D_P2_CLEAN_RT_TRIAGE.tsv outputs/STAGE05D4R5D_SUMMARY.json outputs/GATE_STAGE05D4R5D.json > outputs/STAGE05D4R5D_SHA256.tsv
cat outputs/STAGE05D4R5D_SUMMARY.json
