#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
P=${1:-}
if [[ -z "$P" ]]; then echo "Usage: $0 /absolute/path/to/ART_MathDiscovery_v1" >&2; exit 2; fi
python "$HERE/scripts/resolve_truth.py" --project-root "$P" | tee "$HERE/logs/run_01_truth_resolution.log"
