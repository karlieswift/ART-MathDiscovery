#!/usr/bin/env bash
set -euo pipefail
R5=${1:?R5 package dir}
R4=${2:?R4 package dir}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$R5"
[ -f raw/P1_FULL3D_AFSP.m8 ] || { echo "[FAIL] missing raw/P1_FULL3D_AFSP.m8"; exit 2; }
rows=$(wc -l < raw/P1_FULL3D_AFSP.m8)
echo "[PREFLIGHT] existing_v107_m8_rows=$rows"
echo "[PREFLIGHT] raw m8 query IDs:"
cut -f1 raw/P1_FULL3D_AFSP.m8 | sort | uniq -c || true
python "$HERE/resummarize_r5_v108.py" --r5 "$R5" --r4 "$R4"
echo "[RESULT] query-id map"
cat outputs/STAGE05D4R5_M8_QUERY_ID_MAP.tsv
echo "[RESULT] non-empty full3D hits"
awk -F '\t' 'NR==1 || $10+0>0' outputs/STAGE05D4R5_P1_FULL3D_TRIAGE.tsv
