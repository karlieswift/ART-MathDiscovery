#!/usr/bin/env bash
set -euo pipefail
P=${P:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
B4A=${B4A:-$P/ART_MATH_DISCOVERY_STAGE05F3B4A_SURVIVOR_NATIVE_LOCUS_RECONSTRUCTION_v1_0_0_20261011}
HERE=$(cd "$(dirname "$0")" && pwd)
MODE=${1:-preflight}
WORKERS=${2:-2}
mkdir -p "$HERE/outputs" "$HERE/logs" "$HERE/external/logan_unitigs"
echo '[1/3] context completeness audit'
python "$HERE/scripts/audit_context_completeness.py" --b4a "$B4A" --out "$HERE/outputs"
cat "$HERE/outputs/STAGE05F3B4A2_CONTEXT_AUDIT_GATE.json"
echo '[2/3] Logan unitig availability'
"$HERE/scripts/check_logan_unitigs.sh" "$HERE/outputs/F3B4A2_SOURCE_ACCESSIONS_12.txt" "$HERE/outputs/F3B4A2_LOGAN_UNITIG_AVAILABILITY.tsv"
if [[ "$MODE" == "preflight" ]]; then
  echo 'PASS_STAGE05F3B4A2_CONTEXT_AUDIT_AND_UNITIG_PREFLIGHT_COMPLETE'
  exit 0
fi
if [[ "$MODE" != "download_scan" ]]; then echo "[FATAL] mode must be preflight or download_scan" >&2; exit 2; fi
echo '[3/3] download all 12 frozen-source Logan unitig files and scan exact boundary 31-mers'
"$HERE/scripts/download_logan_unitigs.sh" "$HERE/outputs/F3B4A2_SOURCE_ACCESSIONS_12.txt" "$HERE/external/logan_unitigs" "$WORKERS"
command -v zstdcat >/dev/null 2>&1 || { echo '[FATAL] zstdcat missing' >&2; exit 2; }
python "$HERE/scripts/scan_boundary_anchors.py" --anchors "$HERE/outputs/F3B4A2_BOUNDARY_ANCHORS_31MER.tsv" --unitig-dir "$HERE/external/logan_unitigs" --out "$HERE/outputs/F3B4A2_LOGAN_UNITIG_BOUNDARY_ANCHOR_SCAN.tsv"
sha256sum "$HERE/outputs/"* > "$HERE/outputs/STAGE05F3B4A2_SHA256.tsv"
echo 'PASS_STAGE05F3B4A2_UNITIG_BOUNDARY_ANCHOR_SCAN_COMPLETE'
