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

REUSE_ARGS=()
V100="$P/ART_MATH_DISCOVERY_STAGE05F2E_WHOLEGENOME_ANNOTATION_INDEPENDENT_FAMILY_EXPANSION_v1_0_0_20261009/work"
if [ -s "$V100/STAGE05F2E_WHOLE919_SIXFRAME_INSTANCES.faa" ] && [ -s "$V100/STAGE05F2E_WHOLE919_SIXFRAME_INSTANCES.tsv" ]; then
  echo "[REUSE] completed v1.0.0 whole-genome six-frame materialization: $V100"
  REUSE_ARGS=(--reuse-v100-work "$V100")
else
  echo "[INFO] no reusable v1.0.0 materialization found; rebuilding whole-genome six-frame universe"
fi

python scripts/stage05f2e.py \
  --f2a "$F2A" \
  --f2d "$F2D" \
  --project-root "$P" \
  --threads "$THREADS" \
  --outdir outputs \
  --rawdir raw \
  --work work \
  --logs logs \
  "${REUSE_ARGS[@]}" | tee logs/stage05f2e.log
cat outputs/STAGE05F2E_GATE.json
