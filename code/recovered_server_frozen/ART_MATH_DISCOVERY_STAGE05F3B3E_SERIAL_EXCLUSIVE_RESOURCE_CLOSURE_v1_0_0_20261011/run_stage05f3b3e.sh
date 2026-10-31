#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
GPU_ARG=${2:-auto}
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
D="$P/ART_MATH_DISCOVERY_STAGE05F3B3D_EXACT_PROTOCOL_TECHNICAL_RETRY_v1_0_0_20261011"
V102="$P/ART_MATH_DISCOVERY_STAGE05F3B3C_BOLTZ2_FULL3D_AND_FROZEN_R5C_v1_0_2_20261011"
R5C="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007"
LOCK="$S/frozen/F3B3E_IMPLEMENTATION_LOCK.json"
TARGETS="$S/frozen/F3B3E_RESOURCE_CLOSURE_TARGETS.tsv"
mkdir -p "$S"/{frozen,outputs,logs,serial_boltz_out,raw/foldseek_serial,tmp/foldseek_serial,state}

for x in python boltz conda sha256sum nvidia-smi; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
python - <<'PYENV'
import torch
print('[RUNTIME] torch =',torch.__version__)
print('[RUNTIME] cuda_available =',torch.cuda.is_available())
print('[RUNTIME] cuda_device_count =',torch.cuda.device_count())
if not torch.cuda.is_available(): raise SystemExit('[FATAL] CUDA unavailable')
PYENV
conda run -n art_novelty bash -c 'command -v foldseek >/dev/null' >/dev/null 2>&1 || { echo '[FATAL] foldseek unavailable in art_novelty' >&2; exit 2; }

DG="$D/outputs/STAGE05F3B3D_GATE.json"
DS="$D/outputs/F3B3D_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv"
QM="$V102/frozen/F3B3C_QUERY_MANIFEST.tsv"
RTDB="$R5C/references/clean_rt_afsp"
for x in "$DG" "$DS" "$QM" "$S/scripts/run_serial_exclusive_retry.py" "$S/scripts/summarize_serial_closure.py"; do [[ -s "$x" ]] || { echo "[FATAL] missing $x" >&2; exit 2; }; done
ls "${RTDB}"* >/dev/null 2>&1 || { echo '[FATAL] missing R5C clean-RT DB' >&2; exit 2; }

if [[ ! -s "$LOCK" ]]; then
  echo '[1/7] freeze all-7 serial-exclusive resource closure BEFORE new retry results'
  python - "$DG" "$DS" "$QM" "$S/scripts/run_serial_exclusive_retry.py" "$S/scripts/summarize_serial_closure.py" "$S/run_stage05f3b3e.sh" "$TARGETS" "$LOCK" <<'PY'
import csv,hashlib,json,sys
from pathlib import Path
dg,ds,qm,p1,p2,runner,targets,out=sys.argv[1:9]
def h(p):
 x=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): x.update(b)
 return x.hexdigest()
g=json.load(open(dg))
assert g['gate_pass'] is True
assert g['unresolved_boltz_after_exact_retry_n']==7
assert g['final_r5c_strong_n']==9 and g['final_r5c_suggestive_n']==7 and g['final_r5c_weak_n']==9
with open(ds,newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==32
u=[r for r in rows if r['final_best_tier']=='BOLTZ_TECHNICAL_FAILURE_AFTER_EXACT_RETRY']
assert len(u)==7, len(u)
cols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length']
with open(targets,'w',newline='') as f:
 w=csv.DictWriter(f,delimiter='\t',fieldnames=cols,lineterminator='\n');w.writeheader()
 for r in sorted(u,key=lambda z:int(z['shortlist_rank'])): w.writerow({k:r[k] for k in cols})
lock={
 'stage':'Stage05F3B3E-implementation-freeze','gate':'PASS_F3B3E_SERIAL_EXCLUSIVE_RESOURCE_CLOSURE_FROZEN_PRE_RETRY','gate_pass':True,
 'role':'post-hoc resource-only technical closure of all and only the 7 F3B3D unresolved Boltz CUDA-OOM cases',
 'f3b3d_gate_sha256':h(dg),'f3b3d_final_summary_sha256':h(ds),'query_manifest_sha256':h(qm),'target_tsv_sha256':h(targets),
 'retry_target_n':7,'retry_target_rule':'all and only F3B3D rows with final_best_tier=BOLTZ_TECHNICAL_FAILURE_AFTER_EXACT_RETRY; no length/family/annotation filtering',
 'attempts_per_target_in_f3b3e':1,
 'resource_execution_change':'4 concurrent independent single-GPU jobs -> 1 selected idle physical GPU, strictly serial, fresh Boltz subprocess for each query',
 'boltz_model':'boltz2','msa_mode':'server-assisted','seed':20261007,'recycling_steps':3,'sampling_steps':200,'diffusion_samples':1,'no_kernels':True,'devices_per_job':1,'accelerator':'gpu','output_format':'mmcif',
 'no_single_sequence_fallback':True,'no_length_truncation':True,'no_sequence_or_msa_content_change':True,
 'foldseek_search_evalue':1.0,'foldseek_max_seqs':100,'foldseek_exact_tmscore':1,
 'strong_rule':'alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50','suggestive_rule':'alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30',
 'shortlist_membership_changed':False,'scientific_threshold_changed':False,'primary_or_f3b3d_results_overwritten':False,
 'retry_script_sha256':h(p1),'summarizer_sha256':h(p2),'runner_sha256':h(runner),
 'claim_boundary':'Resource scheduling rescue only. Any newly resolved structure is a technical sensitivity result layered on preserved F3B3C/F3B3D evidence; unresolved cases remain non-negative technical failures.'
}
with open(out,'w') as f: json.dump(lock,f,indent=2);f.write('\n')
print(json.dumps(lock,indent=2))
PY
  find "$R5C/references" -maxdepth 1 -type f -name 'clean_rt_afsp*' -print0 | sort -z | xargs -0 sha256sum > "$S/frozen/R5C_CLEAN_RT_AFSP_SHA256.tsv"
  sha256sum "$LOCK" "$TARGETS" "$S/run_stage05f3b3e.sh" "$S/scripts/run_serial_exclusive_retry.py" "$S/scripts/summarize_serial_closure.py" > "$S/frozen/F3B3E_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$TARGETS" "$S/run_stage05f3b3e.sh" "$S/scripts/run_serial_exclusive_retry.py" "$S/scripts/summarize_serial_closure.py"
else
  echo '[1/7] [SKIP] implementation lock exists'
fi

# Select one idle GPU without changing scientific model settings.
choose_gpu() {
  local arg="$1"
  if [[ "$arg" != auto ]]; then echo "$arg"; return; fi
  local best='' bestfree=-1
  while IFS=',' read -r idx free; do
    idx=$(echo "$idx" | tr -d ' '); free=$(echo "$free" | tr -d ' ')
    [[ "$idx" =~ ^[0-9]+$ && "$free" =~ ^[0-9]+$ ]] || continue
    procs=$(nvidia-smi -i "$idx" --query-compute-apps=pid --format=csv,noheader,nounits 2>/dev/null | grep -E '^[[:space:]]*[0-9]+' || true)
    if [[ -z "$procs" && "$free" -gt "$bestfree" ]]; then best="$idx"; bestfree="$free"; fi
  done < <(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits)
  [[ -n "$best" ]] || { echo '[FATAL] no idle GPU found; do not run resource-closure concurrently with other GPU jobs' >&2; exit 2; }
  echo "$best"
}
GPU_ID=$(choose_gpu "$GPU_ARG")
[[ "$GPU_ID" =~ ^[0-9]+$ ]] || { echo '[FATAL] invalid GPU ID' >&2; exit 2; }
ACTIVE=$(nvidia-smi -i "$GPU_ID" --query-compute-apps=pid --format=csv,noheader,nounits 2>/dev/null | grep -E '^[[:space:]]*[0-9]+' || true)
[[ -z "$ACTIVE" ]] || { echo "[FATAL] GPU $GPU_ID has active compute processes; choose an idle GPU" >&2; exit 2; }
FREE_MB=$(nvidia-smi -i "$GPU_ID" --query-gpu=memory.free --format=csv,noheader,nounits | tr -d ' ')
echo "[RESOURCE] selected_idle_gpu=$GPU_ID free_memory_mb=$FREE_MB mode=SERIAL_EXCLUSIVE"
python - "$GPU_ID" "$FREE_MB" "$S/outputs/F3B3E_RESOURCE_PREFLIGHT.json" <<'PY'
import json,sys
x={'selected_physical_gpu_id':int(sys.argv[1]),'free_memory_mb_at_preflight':int(sys.argv[2]),'active_compute_processes_at_preflight':0,'execution_mode':'SERIAL_EXCLUSIVE_SINGLE_GPU','scientific_model_settings_changed':False}
open(sys.argv[3],'w').write(json.dumps(x,indent=2)+'\n')
PY

echo '[2/7] serial-exclusive Boltz attempt for all 7 unresolved cases'
export BOLTZ_CACHE=${BOLTZ_CACHE:-"$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache"}
mkdir -p "$BOLTZ_CACHE"
python "$S/scripts/run_serial_exclusive_retry.py" "$QM" "$TARGETS" "$S/serial_boltz_out" "$S/logs" "$S/state" "$S/outputs/F3B3E_SERIAL_BOLTZ_STATUS.tsv" "$GPU_ID" | tee "$S/logs/serial_boltz_summary.log"

echo '[3/7] frozen R5C Foldseek for newly resolved structures only'
FSTAT="$S/outputs/F3B3E_SERIAL_FOLDSEEK_STATUS.tsv"
printf 'candidate_id\tfoldseek_status\treturn_code\tm8_path\n' > "$FSTAT"
tail -n +2 "$S/outputs/F3B3E_SERIAL_BOLTZ_STATUS.tsv" | while IFS=$'\t' read -r rank fam cid acc len gpu bst rc elapsed cat cif; do
  [[ "$bst" == PASS* ]] || continue
  out="$S/raw/foldseek_serial/${cid}.m8"; tmp="$S/tmp/foldseek_serial/${cid}"; log="$S/logs/serial_foldseek_${cid}.log"
  rm -rf "$tmp"; rm -f "$out"
  set +e
  conda run -n art_novelty foldseek easy-search "$cif" "$RTDB" "$out" "$tmp" -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 --format-output 'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore' > "$log" 2>&1
  frc=$?
  set -e
  if [[ "$frc" -eq 0 ]]; then fs=PASS; echo "[PASS] serial foldseek $cid"; else fs=TECHNICAL_FAILURE; echo "[WARN] serial foldseek technical failure $cid rc=$frc"; fi
  printf '%s\t%s\t%s\t%s\n' "$cid" "$fs" "$frc" "$out" >> "$FSTAT"
done

echo '[4/7] combine preserved F3B3D evidence with resource-only closure results'
python "$S/scripts/summarize_serial_closure.py" "$DS" "$S/outputs/F3B3E_SERIAL_BOLTZ_STATUS.tsv" "$FSTAT" "$S/raw/foldseek_serial" "$LOCK" "$S/outputs" | tee "$S/logs/summarize_serial_closure.log"

echo '[5/7] unresolved resource-limited cases'
python - "$S/outputs/F3B3E_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv" <<'PY'
import csv,sys
with open(sys.argv[1],newline='') as f:r=list(csv.DictReader(f,delimiter='\t'))
u=[x for x in r if 'UNRESOLVED' in x['f3b3e_final_best_tier'] or 'TECHNICAL_FAILURE' in x['f3b3e_final_best_tier']]
print('resource_or_technical_unresolved_n =',len(u))
for x in u: print(x['shortlist_rank'],x['candidate_id'],x['aa_length'],x['f3b3e_final_best_tier'],x['serial_failure_category'],sep='\t')
PY

echo '[6/7] hashes'
sha256sum "$LOCK" "$TARGETS" "$S/outputs/F3B3E_RESOURCE_PREFLIGHT.json" "$S/outputs/F3B3E_SERIAL_BOLTZ_STATUS.tsv" "$FSTAT" "$S/outputs/F3B3E_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv" "$S/outputs/STAGE05F3B3E_GATE.json" > "$S/outputs/STAGE05F3B3E_SHA256.tsv"

echo '[7/7] compact result'
cat "$S/outputs/STAGE05F3B3E_GATE.json"
echo
echo '=== SERIAL RESOURCE-CLOSURE STATUS 7 ==='
column -t -s $'\t' "$S/outputs/F3B3E_SERIAL_BOLTZ_STATUS.tsv" || true
echo
echo '=== FINAL STRUCTURE/R5C SUMMARY 32 ==='
column -t -s $'\t' "$S/outputs/F3B3E_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv" | head -33 || true
echo
echo 'PASS_STAGE05F3B3E_SERIAL_EXCLUSIVE_RESOURCE_CLOSURE_RUN_COMPLETE'
