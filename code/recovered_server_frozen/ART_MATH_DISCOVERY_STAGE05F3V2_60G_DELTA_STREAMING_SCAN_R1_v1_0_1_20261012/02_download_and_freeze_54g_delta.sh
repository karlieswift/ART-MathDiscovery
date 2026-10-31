#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_MATHDISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
RESUME_FROM=${F3V2_RESUME_FROM:-$P/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_v1_0_0_20261011}
for STR in A B C D; do
  [[ -s "$S/inputs/F3V2_DELTA_CANDIDATES_${STR}.tsv" ]] || { echo '[FATAL] run 01_select_logan_delta_candidates.sh first' >&2; exit 2; }
done
[[ -s "$S/inputs/OLD_F3_ACCESSIONS.txt" ]] || { echo '[FATAL] missing OLD_F3_ACCESSIONS.txt' >&2; exit 2; }
mkdir -p "$S"/{data/contigs,work/download,logs,frozen}
python "$S/scripts/acquire_54g_r1.py" --stage "$S" --resume-from "$RESUME_FROM"
