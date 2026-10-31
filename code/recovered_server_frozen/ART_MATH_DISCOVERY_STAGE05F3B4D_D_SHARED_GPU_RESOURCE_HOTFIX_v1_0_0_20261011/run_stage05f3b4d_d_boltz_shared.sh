#!/usr/bin/env bash
set -euo pipefail

GPU_ARG=${1:-auto}
MIN_FREE_MB=${F3B4D_D_MIN_FREE_MB:-60000}
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
BASE="$P/ART_MATH_DISCOVERY_STAGE05F3B4D_SEVEN_FAMILY_A_TO_E_CLOSEOUT_v1_0_0_20261011"
HOT=$(cd "$(dirname "$0")" && pwd)
TARGETS="$BASE/outputs/F3B4D_D_STRUCTURE_EXTENSION_TARGETS_6.tsv"
STATUS="$BASE/outputs/F3B4D_D_BOLTZ_STATUS.tsv"
RUNNER="$BASE/scripts/run_d_boltz.py"
LOCK="$HOT/frozen/F3B4D_D_SHARED_GPU_RESOURCE_HOTFIX_LOCK.json"

for x in python boltz nvidia-smi sha256sum; do
  command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }
done
[[ -s "$TARGETS" ]] || { echo '[FATAL] base Stage05F3B4D prepare outputs missing' >&2; exit 2; }
[[ -s "$RUNNER" ]] || { echo '[FATAL] base D Boltz runner missing' >&2; exit 2; }

# Result-blind hotfix only. The earlier idle-GPU attempt aborted before any candidate run.
if [[ -s "$STATUS" ]]; then
  echo "[FATAL] $STATUS already exists; shared-GPU hotfix is only authorized pre-D-result" >&2
  exit 2
fi

EXPECTED_RUNNER_SHA=616ca06b2d2fff7f7e709d98673122e44f217c5aacd437147768caf740c757ef
ACTUAL_RUNNER_SHA=$(sha256sum "$RUNNER" | awk '{print $1}')
[[ "$ACTUAL_RUNNER_SHA" == "$EXPECTED_RUNNER_SHA" ]] || {
  echo "[FATAL] base run_d_boltz.py hash mismatch: $ACTUAL_RUNNER_SHA" >&2
  exit 2
}

python - "$TARGETS" <<'PY'
import csv,sys
r=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
assert len(r)==6, len(r)
assert [x['family_rep_id'] for x in r]==[
 'F3BI37356863','F3BI37262050','F3BI27373723',
 'F3BI55368798','F3BI05882631','F3BI52178020']
print('[PASS] frozen six-target manifest verified')
PY

python - <<'PY'
import torch
print('[RUNTIME] torch =',torch.__version__)
print('[RUNTIME] cuda_available =',torch.cuda.is_available())
print('[RUNTIME] cuda_device_count =',torch.cuda.device_count())
if not torch.cuda.is_available(): raise SystemExit('[FATAL] CUDA unavailable')
PY

choose_gpu() {
  local arg="$1"
  if [[ "$arg" != auto ]]; then
    echo "$arg"
    return
  fi
  local best='' bestfree=-1
  while IFS=',' read -r idx free; do
    idx=$(echo "$idx" | tr -d ' ')
    free=$(echo "$free" | tr -d ' ')
    [[ "$idx" =~ ^[0-9]+$ && "$free" =~ ^[0-9]+$ ]] || continue
    if (( free > bestfree )); then best="$idx"; bestfree="$free"; fi
  done < <(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits)
  [[ -n "$best" ]] || { echo '[FATAL] cannot resolve GPU' >&2; exit 2; }
  echo "$best"
}

GPU=$(choose_gpu "$GPU_ARG")
FREE=$(nvidia-smi -i "$GPU" --query-gpu=memory.free --format=csv,noheader,nounits | tr -d ' ')
TOTAL=$(nvidia-smi -i "$GPU" --query-gpu=memory.total --format=csv,noheader,nounits | tr -d ' ')
UTIL=$(nvidia-smi -i "$GPU" --query-gpu=utilization.gpu --format=csv,noheader,nounits | tr -d ' %')
[[ "$FREE" =~ ^[0-9]+$ ]] || { echo '[FATAL] invalid free-memory query' >&2; exit 2; }
if (( FREE < MIN_FREE_MB )); then
  echo "[FATAL] GPU $GPU free_memory_mb=$FREE < frozen shared-GPU launch floor $MIN_FREE_MB" >&2
  exit 2
fi
ACTIVE=$(nvidia-smi -i "$GPU" --query-compute-apps=pid,used_memory --format=csv,noheader,nounits 2>/dev/null | sed '/^[[:space:]]*$/d' || true)
ACTIVE_ONE_LINE=$(printf '%s' "$ACTIVE" | tr '\n' ';' | sed 's/;*$//')

echo "[RESOURCE] selected_gpu=$GPU free_memory_mb=$FREE total_memory_mb=$TOTAL gpu_util_pct=$UTIL mode=SERIAL_SHARED_GPU"
echo "[RESOURCE] preexisting_compute_apps=${ACTIVE_ONE_LINE:-NONE}"

mkdir -p "$HOT/frozen" "$HOT/logs" "$BASE/structure_boltz_out" "$BASE/logs" "$BASE/state"

if [[ ! -s "$LOCK" ]]; then
  python - "$TARGETS" "$RUNNER" "$LOCK" "$GPU" "$FREE" "$TOTAL" "$UTIL" "$MIN_FREE_MB" "$ACTIVE_ONE_LINE" <<'PY'
import csv,hashlib,json,sys,datetime
T,R,O,gpu,free,total,util,minfree,active=sys.argv[1:10]
def h(p):
    z=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): z.update(b)
    return z.hexdigest()
rows=list(csv.DictReader(open(T),delimiter='\t'))
obj={
 'stage':'Stage05F3B4D-D-resource-hotfix',
 'gate':'PASS_F3B4D_D_SHARED_GPU_RESOURCE_POLICY_FROZEN_PRE_RESULTS',
 'gate_pass':True,
 'hotfix_role':'RESULT_BLIND_RESOURCE_SCHEDULING_ONLY',
 'reason':'idle-exclusive launch rule was unnecessarily conservative for six short 115-434 aa targets; previous D attempt aborted before any candidate result',
 'target_n':len(rows),
 'target_manifest_sha256':h(T),
 'base_boltz_runner_sha256':h(R),
 'boltz_scientific_settings_changed':False,
 'target_membership_changed':False,
 'sequence_changed':False,
 'structure_threshold_changed':False,
 'resource_policy_old':'SERIAL_EXCLUSIVE_IDLE_GPU_ONLY',
 'resource_policy_new':'SERIAL_SHARED_GPU_ALLOWED_IF_LAUNCH_FREE_MEMORY_FLOOR_IS_MET',
 'shared_gpu_launch_min_free_memory_mb':int(minfree),
 'selected_physical_gpu_id':int(gpu),
 'launch_free_memory_mb':int(free),
 'launch_total_memory_mb':int(total),
 'launch_gpu_utilization_pct':int(util),
 'preexisting_compute_apps':active or 'NONE',
 'failure_handling':'Any CUDA OOM or other failure remains technical/resource-censored and may be retried later without changing scientific settings.',
 'claim_boundary':'Resource scheduling hotfix only; sharing a busy GPU changes runtime/technical-failure risk, not the frozen Boltz/R5C scientific protocol.'
}
open(O,'w').write(json.dumps(obj,indent=2)+'\n')
print(json.dumps(obj,indent=2))
PY
else
  echo '[SKIP] shared-GPU resource hotfix lock already exists'
fi

export BOLTZ_CACHE=${BOLTZ_CACHE:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache}
mkdir -p "$BOLTZ_CACHE"
python "$RUNNER" \
  "$TARGETS" \
  "$BASE/structure_boltz_out" \
  "$BASE/logs" \
  "$BASE/state" \
  "$STATUS" \
  "$GPU" \
  | tee "$HOT/logs/d_shared_boltz_summary.log"

echo PASS_STAGE05F3B4D_D_SHARED_GPU_BOLTZ_RUN_COMPLETE
