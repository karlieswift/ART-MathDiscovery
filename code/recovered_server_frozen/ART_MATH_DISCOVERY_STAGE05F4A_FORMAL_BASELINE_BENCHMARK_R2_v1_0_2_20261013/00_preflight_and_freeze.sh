#!/usr/bin/env bash
set -euo pipefail
P=${1:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/logs" "$HERE/frozen" "$HERE/outputs" "$HERE/work"
python "$HERE/scripts/inventory_and_freeze.py" "$P" | tee "$HERE/logs/run_00_preflight.log"
python "$HERE/scripts/suggest_truth_sources.py" "$P" | tee "$HERE/logs/run_00_truth_source_inventory.log"
