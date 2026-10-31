#!/usr/bin/env bash
set -euo pipefail
# usage: ./run_stage05f2e.sh F2A_PATH F2D_PATH PROJECT_ROOT [THREADS]
F2A=${1:?F2A path required}
F2D=${2:?F2D path required}
P=${3:?project root required}
THREADS=${4:-32}
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
rm -rf outputs raw work logs
mkdir -p outputs raw work logs
python scripts/stage05f2e.py --selftest
python scripts/stage05f2e.py \
  --f2a "$F2A" \
  --f2d "$F2D" \
  --project-root "$P" \
  --threads "$THREADS" \
  --outdir outputs \
  --rawdir raw \
  --work work \
  --logs logs | tee logs/stage05f2e.log
cat outputs/STAGE05F2E_GATE.json
