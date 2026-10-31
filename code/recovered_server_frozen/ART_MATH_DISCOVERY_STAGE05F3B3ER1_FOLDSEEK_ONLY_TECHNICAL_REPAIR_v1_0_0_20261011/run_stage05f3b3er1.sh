#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
E="$P/ART_MATH_DISCOVERY_STAGE05F3B3E_SERIAL_EXCLUSIVE_RESOURCE_CLOSURE_v1_0_0_20261011"
R5C="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007"
RTDB="$R5C/references/clean_rt_afsp"
LOCK="$S/frozen/F3B3ER1_IMPLEMENTATION_LOCK.json"
TARGETS="$S/frozen/F3B3ER1_FOLDSEEK_REPAIR_TARGETS.tsv"
ESTAT="$E/outputs/F3B3E_SERIAL_BOLTZ_STATUS.tsv"
ESUM="$E/outputs/F3B3E_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv"
EGATE="$E/outputs/STAGE05F3B3E_GATE.json"
mkdir -p "$S"/{outputs,logs,frozen,raw/foldseek_r1,tmp/foldseek_r1}

for x in python sha256sum foldseek; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
python - <<'PY'
import os
print('[RUNTIME] CONDA_DEFAULT_ENV =',os.environ.get('CONDA_DEFAULT_ENV',''))
print('[RUNTIME] CONDA_PREFIX =',os.environ.get('CONDA_PREFIX',''))
PY
if [[ "${CONDA_DEFAULT_ENV:-}" != *art_novelty* && "${CONDA_PREFIX:-}" != *art_novelty* ]]; then
  echo '[FATAL] activate art_novelty before running; this stage deliberately avoids conda-run wrapper' >&2
  exit 2
fi
for x in "$ESTAT" "$ESUM" "$EGATE" "$S/scripts/summarize_foldseek_repair.py"; do [[ -s "$x" ]] || { echo "[FATAL] missing $x" >&2; exit 2; }; done
ls "${RTDB}"* >/dev/null 2>&1 || { echo '[FATAL] missing frozen R5C DB' >&2; exit 2; }

if [[ ! -s "$LOCK" ]]; then
  echo '[1/5] freeze Foldseek-only technical repair BEFORE repaired Foldseek results'
  python - "$ESTAT" "$ESUM" "$EGATE" "$TARGETS" "$LOCK" "$S/run_stage05f3b3er1.sh" "$S/scripts/summarize_foldseek_repair.py" <<'PY'
import csv,hashlib,json,sys
estat,esum,egate,targets,lock,runner,summ=sys.argv[1:8]
def h(p):
 x=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): x.update(b)
 return x.hexdigest()
g=json.load(open(egate)); assert g['gate_pass'] is True
assert g['serial_exclusive_boltz_success_n']==5 and g['foldseek_technical_failure_n']==5
with open(esum,newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
t=[r for r in rows if r['f3b3e_final_best_tier']=='FOLDSEEK_TECHNICAL_FAILURE_AFTER_SERIAL_RESOURCE_CLOSURE']
assert len(t)==5,len(t)
cols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','serial_cif_path']
with open(targets,'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n'); w.writeheader()
 for r in sorted(t,key=lambda z:int(z['shortlist_rank'])):
  assert r['serial_exclusive_boltz_status'].startswith('PASS')
  w.writerow({k:r[k] for k in cols})
lockd={
 'stage':'Stage05F3B3E-R1-implementation-freeze',
 'gate':'PASS_F3B3ER1_FOLDSEEK_ONLY_REPAIR_FROZEN_PRE_RESULTS','gate_pass':True,
 'role':'post-hoc technical Foldseek handoff/runtime repair only; no Boltz recomputation',
 'target_rule':'all and only the five F3B3E rows with PASS_SERIAL_EXCLUSIVE Boltz and FOLDSEEK_TECHNICAL_FAILURE_AFTER_SERIAL_RESOURCE_CLOSURE',
 'target_n':5,'boltz_recomputed':False,'cif_reused_exactly':True,
 'foldseek_runtime':'direct foldseek executable from explicitly activated art_novelty environment; no conda-run wrapper',
 'foldseek_evalue':1.0,'foldseek_max_seqs':100,'foldseek_threads_cli_argument':'user supplied; default 32','foldseek_exact_tmscore':1,
 'format_output':'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore',
 'strong_rule':'alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50',
 'suggestive_rule':'alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30',
 'scientific_threshold_changed':False,'shortlist_changed':False,'sequence_changed':False,'structure_changed':False,
 'f3b3e_gate_sha256':h(egate),'f3b3e_summary_sha256':h(esum),'f3b3e_serial_status_sha256':h(estat),
 'targets_sha256':h(targets),'runner_sha256':h(runner),'summarizer_sha256':h(summ),
 'claim_boundary':'Technical Foldseek repair only. Repaired STRONG/SUGGESTIVE remains frozen clean-RT structural compatibility, not biochemical RT/ART identity or novelty.'}
with open(lock,'w') as f: json.dump(lockd,f,indent=2);f.write('\n')
print(json.dumps(lockd,indent=2))
PY
  sha256sum "$LOCK" "$TARGETS" > "$S/frozen/F3B3ER1_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$TARGETS" "$S/run_stage05f3b3er1.sh" "$S/scripts/summarize_foldseek_repair.py"
else
  echo '[1/5] [SKIP] implementation lock exists'
fi

echo '[2/5] direct frozen R5C Foldseek on the 5 existing serial-exclusive CIFs only'
FSTAT="$S/outputs/F3B3ER1_FOLDSEEK_STATUS.tsv"
printf 'candidate_id\tfoldseek_status\treturn_code\tm8_path\tlog_path\n' > "$FSTAT"
tail -n +2 "$TARGETS" | while IFS=$'\t' read -r rank fam cid acc len cif; do
  [[ -s "$cif" ]] || { echo "[FATAL] missing CIF $cid: $cif" >&2; exit 2; }
  out="$S/raw/foldseek_r1/${cid}.m8"; tmp="$S/tmp/foldseek_r1/${cid}"; log="$S/logs/foldseek_r1_${cid}.log"
  rm -rf "$tmp"; rm -f "$out"
  set +e
  foldseek easy-search "$cif" "$RTDB" "$out" "$tmp" -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 --format-output 'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore' > "$log" 2>&1
  rc=$?
  set -e
  if [[ "$rc" -eq 0 && -e "$out" ]]; then
    fs=PASS; echo "[PASS] foldseek repair $cid"
  else
    fs=TECHNICAL_FAILURE; echo "[WARN] foldseek repair failed $cid rc=$rc"
    tail -20 "$log" || true
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$cid" "$fs" "$rc" "$out" "$log" >> "$FSTAT"
done

echo '[3/5] integrate repaired Foldseek results without modifying F3B3E outputs'
python "$S/scripts/summarize_foldseek_repair.py" "$ESUM" "$FSTAT" "$S/raw/foldseek_r1" "$LOCK" "$S/outputs"

echo '[4/5] hashes'
sha256sum "$LOCK" "$TARGETS" "$FSTAT" "$S/outputs/F3B3ER1_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv" "$S/outputs/STAGE05F3B3ER1_GATE.json" > "$S/outputs/STAGE05F3B3ER1_SHA256.tsv"

echo '[5/5] compact result'
cat "$S/outputs/STAGE05F3B3ER1_GATE.json"
echo
echo '=== FOLDSEEK REPAIR STATUS 5 ==='
column -t -s $'\t' "$FSTAT" || true
echo
echo '=== REPAIRED TARGET TIERS ==='
python - "$S/outputs/F3B3ER1_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv" <<'PY'
import csv,sys
r=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
for x in r:
 if x.get('b3er1_resolution_source')=='FOLDSEEK_ONLY_TECHNICAL_REPAIR':
  print(x['shortlist_rank'],x['candidate_id'],x['aa_length'],x['b3er1_final_best_tier'],x.get('best_alntmscore',''),x.get('best_qcov',''),x.get('best_tcov',''),sep='\t')
PY

echo PASS_STAGE05F3B3ER1_FOLDSEEK_ONLY_TECHNICAL_REPAIR_RUN_COMPLETE
