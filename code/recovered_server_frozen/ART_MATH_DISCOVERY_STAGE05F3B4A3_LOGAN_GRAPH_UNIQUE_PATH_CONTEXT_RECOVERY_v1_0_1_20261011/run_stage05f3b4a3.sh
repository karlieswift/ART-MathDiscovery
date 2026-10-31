#!/usr/bin/env bash
set -euo pipefail
P=${P:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
B4A=${B4A:-$P/ART_MATH_DISCOVERY_STAGE05F3B4A_SURVIVOR_NATIVE_LOCUS_RECONSTRUCTION_v1_0_0_20261011}
B4A2=${B4A2:-$P/ART_MATH_DISCOVERY_STAGE05F3B4A2_CONTEXT_COMPLETENESS_AND_LOGAN_UNITIG_PREFLIGHT_v1_0_0_20261011}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE"/{outputs,logs,raw,work}
command -v zstdcat >/dev/null 2>&1 || { echo '[FATAL] zstdcat missing' >&2; exit 2; }
[[ -f "$B4A2/outputs/F3B4A2_LOGAN_UNITIG_BOUNDARY_ANCHOR_SCAN.tsv" ]] || { echo '[FATAL] missing B4A2 anchor scan' >&2; exit 2; }
[[ -f "$B4A/outputs/F3B4A_SURVIVOR_NATIVE_LOCUS_MAP.tsv" ]] || { echo '[FATAL] missing B4A locus map' >&2; exit 2; }
echo '[0/3] selftest'
"$HERE/tests/run_selftest.sh"
echo '[1/3] freeze unique-path graph rule BEFORE graph extension results'
python "$HERE/scripts/freeze_implementation.py" "$HERE" "$B4A2" "$B4A"
echo '[2/3] build Logan graph metadata accession-by-accession and recover only unbranched source context'
python "$HERE/scripts/graph_unique_path_recovery.py" --b4a "$B4A" --b4a2 "$B4A2" --out "$HERE/outputs" --work "$HERE/work"
echo '[3/3] compact closeout'
cat "$HERE/outputs/STAGE05F3B4A3_GATE.json"
echo 'PASS_STAGE05F3B4A3_RUN_COMPLETE'
