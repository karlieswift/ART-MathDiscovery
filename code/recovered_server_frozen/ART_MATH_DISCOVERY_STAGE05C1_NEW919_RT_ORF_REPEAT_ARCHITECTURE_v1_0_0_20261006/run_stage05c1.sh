#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
A05="${1:?usage: run_stage05c1.sh STAGE05A_ROOT STAGE05B_ROOT [CPUS=4]}"
B05="${2:?usage: run_stage05c1.sh STAGE05A_ROOT STAGE05B_ROOT [CPUS=4]}"
CPUS="${3:-4}"
cd "$ROOT"
mkdir -p logs
"${PYTHON:-python}" -u scripts/stage05c1.py run --source-root "$A05" --score-root "$B05" --cpus "$CPUS" 2>&1 | tee -a logs/run_stage05c1.log
