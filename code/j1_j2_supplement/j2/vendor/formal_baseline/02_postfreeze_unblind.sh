#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
TRUTH=${1:-$HERE/frozen/BASELINE_TRUTH_SEAL.tsv}
python "$HERE/scripts/postfreeze_unblind.py" --truth-tsv "$TRUTH" | tee "$HERE/logs/run_02_postfreeze_unblind.log"
