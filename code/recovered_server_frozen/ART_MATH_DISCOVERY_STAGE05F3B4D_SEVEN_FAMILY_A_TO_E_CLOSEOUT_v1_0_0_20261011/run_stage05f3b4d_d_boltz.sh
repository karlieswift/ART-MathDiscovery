#!/usr/bin/env bash
set -euo pipefail
GPU_ARG=${1:-auto}
S=$(cd "$(dirname "$0")" && pwd)
TARGETS="$S/outputs/F3B4D_D_STRUCTURE_EXTENSION_TARGETS_6.tsv"
[[ -s "$TARGETS" ]] || { echo '[FATAL] run prepare first'; exit 2; }
for x in python boltz nvidia-smi; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x"; exit 2; }; done
python - <<'PY'
import torch
print('[RUNTIME] torch =',torch.__version__)
print('[RUNTIME] cuda_available =',torch.cuda.is_available())
print('[RUNTIME] cuda_device_count =',torch.cuda.device_count())
if not torch.cuda.is_available(): raise SystemExit('[FATAL] CUDA unavailable')
PY
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
  [[ -n "$best" ]] || { echo '[FATAL] no idle GPU; D requires a serial-exclusive idle GPU' >&2; exit 2; }
  echo "$best"
}
GPU=$(choose_gpu "$GPU_ARG")
ACTIVE=$(nvidia-smi -i "$GPU" --query-compute-apps=pid --format=csv,noheader,nounits 2>/dev/null | grep -E '^[[:space:]]*[0-9]+' || true)
[[ -z "$ACTIVE" ]] || { echo "[FATAL] GPU $GPU has active compute processes"; exit 2; }
FREE=$(nvidia-smi -i "$GPU" --query-gpu=memory.free --format=csv,noheader,nounits | tr -d ' ')
echo "[RESOURCE] selected_idle_gpu=$GPU free_memory_mb=$FREE mode=SERIAL_EXCLUSIVE"
export BOLTZ_CACHE=${BOLTZ_CACHE:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache}
mkdir -p "$BOLTZ_CACHE" "$S/structure_boltz_out" "$S/logs" "$S/state"
python "$S/scripts/run_d_boltz.py" "$TARGETS" "$S/structure_boltz_out" "$S/logs" "$S/state" "$S/outputs/F3B4D_D_BOLTZ_STATUS.tsv" "$GPU" | tee "$S/logs/d_boltz_summary.log"
echo PASS_STAGE05F3B4D_D_BOLTZ_RUN_COMPLETE
