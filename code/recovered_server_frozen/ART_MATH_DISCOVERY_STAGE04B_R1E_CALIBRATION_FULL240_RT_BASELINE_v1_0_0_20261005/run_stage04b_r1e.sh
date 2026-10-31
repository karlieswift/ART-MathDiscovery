#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
P="${1:?Usage: run_stage04b_r1e.sh PROJECT_ROOT [CPUS=4]}"
CPUS="${2:-4}"
A04="$P/ART_MATH_DISCOVERY_STAGE04A_L4_BROAD_UNIVERSE_FROZEN_DISCOVERY_v1_0_0_20260928"
D2="$P/ART_MATH_DISCOVERY_STAGE04B_R1D2_RT_ORF_ARRAY_EVIDENCE_v1_0_0_20261005"
mkdir -p logs
python -u scripts/run_stage.py --stage04a "$A04" --r1d2 "$D2" --cpus "$CPUS" 2>&1 | tee logs/run_stage04b_r1e.log
