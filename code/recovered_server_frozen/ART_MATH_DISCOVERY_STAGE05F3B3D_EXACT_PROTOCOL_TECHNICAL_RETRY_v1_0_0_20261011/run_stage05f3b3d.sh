#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
NGPU=${2:-4}
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
V102="$P/ART_MATH_DISCOVERY_STAGE05F3B3C_BOLTZ2_FULL3D_AND_FROZEN_R5C_v1_0_2_20261011"
V103="$P/ART_MATH_DISCOVERY_STAGE05F3B3C_BOLTZ2_FULL3D_AND_FROZEN_R5C_v1_0_3_20261011"
R5C="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007"
LOCK="$S/frozen/F3B3D_IMPLEMENTATION_LOCK.json"

mkdir -p "$S"/{frozen,outputs,logs,retry_boltz_out,raw/foldseek_retry,tmp/foldseek_retry}

for x in python boltz conda sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
python - <<'PYENV'
import torch
print('[RUNTIME] torch =',torch.__version__)
print('[RUNTIME] cuda_available =',torch.cuda.is_available())
print('[RUNTIME] cuda_device_count =',torch.cuda.device_count())
if not torch.cuda.is_available(): raise SystemExit('[FATAL] CUDA unavailable')
PYENV
conda run -n art_novelty bash -c 'command -v foldseek >/dev/null' >/dev/null 2>&1 || { echo '[FATAL] foldseek unavailable in art_novelty' >&2; exit 2; }
[[ "$NGPU" -ge 1 && "$NGPU" -le 4 ]] || { echo '[FATAL] NGPU must be 1..4' >&2; exit 2; }

PG="$V103/outputs/STAGE05F3B3C_GATE.json"
PS="$V103/outputs/F3B3C_STRUCTURE_R5C_SUMMARY_32.tsv"
BS="$V102/outputs/F3B3C_BOLTZ_STATUS.tsv"
QM="$V102/frozen/F3B3C_QUERY_MANIFEST.tsv"
RTDB="$R5C/references/clean_rt_afsp"
for x in "$PG" "$PS" "$BS" "$QM" "$S/scripts/diagnose_primary_failures.py" "$S/scripts/run_exact_retry.py" "$S/scripts/summarize_retry.py"; do [[ -s "$x" ]] || { echo "[FATAL] missing $x" >&2; exit 2; }; done
ls "${RTDB}"* >/dev/null 2>&1 || { echo '[FATAL] missing R5C DB' >&2; exit 2; }

if [[ ! -s "$LOCK" ]]; then
  echo '[1/7] freeze exact-protocol retry rule BEFORE failure-log diagnosis'
  find "$R5C/references" -maxdepth 1 -type f -name 'clean_rt_afsp*' -print0 | sort -z | xargs -0 sha256sum > "$S/frozen/R5C_CLEAN_RT_AFSP_SHA256.tsv"
  python - "$PG" "$PS" "$BS" "$QM" "$S/frozen/R5C_CLEAN_RT_AFSP_SHA256.tsv" "$S/scripts/diagnose_primary_failures.py" "$S/scripts/run_exact_retry.py" "$S/scripts/summarize_retry.py" "$S/run_stage05f3b3d.sh" "$LOCK" <<'PY'
import csv,hashlib,json,sys
pg,ps,bs,qm,dbm,p1,p2,p3,runner,out=sys.argv[1:11]
def h(p):
 x=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): x.update(b)
 return x.hexdigest()
g=json.load(open(pg))
assert g['gate_pass'] is True and g['structure_handoff_family_n']==32
assert g['boltz_success_n']==22 and g['boltz_technical_failure_n']==10
assert g['r5c_strong_or_suggestive_survivor_n']==14
with open(bs,newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
fails=[r['candidate_id'] for r in rows if r['boltz_status'].strip()=='TECHNICAL_FAILURE']
assert len(rows)==32 and len(fails)==10
lock={
 'stage':'Stage05F3B3D-implementation-freeze',
 'gate':'PASS_F3B3D_EXACT_PROTOCOL_TECHNICAL_RETRY_FROZEN_BEFORE_FAILURE_DIAGNOSIS',
 'gate_pass':True,
 'role':'post-primary technical rescue only',
 'primary_gate_sha256':h(pg),'primary_summary_sha256':h(ps),'primary_boltz_status_sha256':h(bs),'query_manifest_sha256':h(qm),'r5c_db_manifest_sha256':h(dbm),
 'primary_boltz_success_n':22,'primary_boltz_technical_failure_n':10,'retry_target_n':10,
 'retry_target_rule':'all and only primary F3B3C Boltz TECHNICAL_FAILURE rows; no length/cause/content filtering',
 'retry_attempts_per_target':1,
 'diagnostic_category_used_for_retry_selection':False,
 'boltz_model':'boltz2','msa_mode':'server-assisted','seed':20261007,'recycling_steps':3,'sampling_steps':200,'diffusion_samples':1,'no_kernels':True,'devices_per_job':1,'accelerator':'gpu','output_format':'mmcif',
 'no_single_sequence_fallback':True,'no_length_truncation':True,
 'foldseek_search_evalue':1.0,'foldseek_max_seqs':100,'foldseek_exact_tmscore':1,
 'strong_rule':'alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50','suggestive_rule':'alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30',
 'shortlist_membership_changed':False,'scientific_threshold_changed':False,'primary_results_overwritten':False,
 'diagnose_script_sha256':h(p1),'retry_script_sha256':h(p2),'summarizer_sha256':h(p3),'runner_sha256':h(runner),
 'claim_boundary':'Exact-protocol technical retry only. Primary F3B3C results remain separately preserved. Retry cannot establish biochemical RT/ART identity or novelty.'
}
with open(out,'w') as f: json.dump(lock,f,indent=2);f.write('\n')
print(json.dumps(lock,indent=2))
PY
  sha256sum "$LOCK" "$S/run_stage05f3b3d.sh" "$S/scripts/diagnose_primary_failures.py" "$S/scripts/run_exact_retry.py" "$S/scripts/summarize_retry.py" > "$S/frozen/F3B3D_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$S/run_stage05f3b3d.sh" "$S/scripts/diagnose_primary_failures.py" "$S/scripts/run_exact_retry.py" "$S/scripts/summarize_retry.py"
else
  echo '[1/7] [SKIP] implementation lock exists'
fi

echo '[2/7] diagnose the 10 primary failures descriptively'
python "$S/scripts/diagnose_primary_failures.py" "$BS" "$V102/logs" "$PS" "$S/outputs" | tee "$S/logs/diagnose_primary_failures.log"

echo '[3/7] exact-protocol Boltz retry of ALL 10 primary technical failures'
export BOLTZ_CACHE=${BOLTZ_CACHE:-"$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache"}
mkdir -p "$BOLTZ_CACHE"
python "$S/scripts/run_exact_retry.py" "$QM" "$BS" "$S/retry_boltz_out" "$S/logs" "$S/outputs/F3B3D_RETRY_BOLTZ_STATUS.tsv" "$NGPU" | tee "$S/logs/retry_boltz_summary.log"

echo '[4/7] frozen R5C Foldseek for retry successes only'
FSTAT="$S/outputs/F3B3D_RETRY_FOLDSEEK_STATUS.tsv"
printf 'candidate_id\tfoldseek_status\treturn_code\tm8_path\n' > "$FSTAT"
tail -n +2 "$S/outputs/F3B3D_RETRY_BOLTZ_STATUS.tsv" | while IFS=$'\t' read -r rank fam cid acc len gpu bst rc elapsed cif; do
  [[ "$bst" == PASS* ]] || continue
  out="$S/raw/foldseek_retry/${cid}.m8"; tmp="$S/tmp/foldseek_retry/${cid}"; log="$S/logs/retry_foldseek_${cid}.log"
  rm -rf "$tmp"; rm -f "$out"
  set +e
  conda run -n art_novelty foldseek easy-search "$cif" "$RTDB" "$out" "$tmp" -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 --format-output 'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore' > "$log" 2>&1
  frc=$?
  set -e
  if [[ "$frc" -eq 0 ]]; then fs=PASS; echo "[PASS] retry foldseek $cid"; else fs=TECHNICAL_FAILURE; echo "[WARN] retry foldseek technical failure $cid rc=$frc"; fi
  printf '%s\t%s\t%s\t%s\n' "$cid" "$fs" "$frc" "$out" >> "$FSTAT"
done

echo '[5/7] combine preserved primary results with separate exact-protocol retry results'
python "$S/scripts/summarize_retry.py" "$PS" "$S/outputs/F3B3D_RETRY_BOLTZ_STATUS.tsv" "$FSTAT" "$S/raw/foldseek_retry" "$LOCK" "$S/outputs" | tee "$S/logs/summarize_retry.log"

echo '[6/7] hashes'
sha256sum "$LOCK" "$S/outputs/F3B3D_PRIMARY_FAILURE_DIAGNOSTICS.tsv" "$S/outputs/F3B3D_RETRY_BOLTZ_STATUS.tsv" "$FSTAT" "$S/outputs/F3B3D_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv" "$S/outputs/STAGE05F3B3D_GATE.json" > "$S/outputs/STAGE05F3B3D_SHA256.tsv"

echo '[7/7] compact result'
cat "$S/outputs/STAGE05F3B3D_GATE.json"
echo
echo '=== PRIMARY FAILURE DIAGNOSIS 10 ==='
column -t -s $'\t' "$S/outputs/F3B3D_PRIMARY_FAILURE_DIAGNOSTICS.tsv" | cut -c1-240 || true
echo
echo '=== FINAL STRUCTURE/R5C SUMMARY 32 ==='
column -t -s $'\t' "$S/outputs/F3B3D_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv" | head -33 || true
echo
echo 'PASS_STAGE05F3B3D_EXACT_PROTOCOL_TECHNICAL_RETRY_COMPLETE'
