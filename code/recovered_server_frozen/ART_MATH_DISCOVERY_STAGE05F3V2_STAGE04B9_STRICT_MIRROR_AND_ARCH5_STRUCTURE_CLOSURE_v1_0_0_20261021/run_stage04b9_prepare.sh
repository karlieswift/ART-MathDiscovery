#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root}; HERE=$(cd "$(dirname "$0")" && pwd)
HIST="$P/ART_MATH_DISCOVERY_STAGE05F3B3A_POSTFREEZE_PRIORITY_AND_CATALYTIC_TOPOLOGY_v1_0_0_20261011"
R3="$P/ART_MATH_DISCOVERY_STAGE05F3B1R3_QUERY_NAMESPACE_RESCUE_AND_F3B2_FAMILY_COMPRESSION_v1_0_0_20261011"
R4="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_v1_0_0_20261011"
mkdir -p "$HERE"/{outputs,work/historical_replay,logs,inputs}
"$HERE/run_selftest.sh"
HISTSCRIPT="$HIST/scripts/stage05f3b3a.py"; EXACT="$HERE/frozen/historical_stage05f3b3a_exact.py"
[[ -s "$HISTSCRIPT" && -s "$HIST/frozen/F3B3A_IMPLEMENTATION_LOCK.json" && -s "$HIST/outputs/F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv" ]] || { echo HOLD historical F3B3A missing; exit 2; }
[[ $(sha256sum "$HISTSCRIPT"|awk '{print $1}') == d222f638c627566909765012874c7aaa531a1a79b5e67a0d73f726f184a30ed8 ]] || { echo HOLD historical F3B3A SHA mismatch; exit 2; }
rm -rf "$HERE/work/historical_replay";mkdir -p "$HERE/work/historical_replay"
python "$EXACT" \
 "$R4/outputs/STAGE05F3B2_R4_GATE.json" \
 "$R4/outputs/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv" \
 "$R4/outputs/F3B2_R4_FAMILY_SUMMARY_WITH_SOURCE.tsv" \
 "$R3/outputs/STAGE05F3B1R3_FORMAL_INSTANCE_IDS.txt" \
 "$R3/raw/F3B1R3_NAMESPACE_INSTANCES.faa" \
 "$R3/outputs/STAGE05F3B1R3_INSTANCE_RVT_PAIR_HITS.tsv.gz" \
 "$R3/raw/F3B1R3_RVT_ALL.domtblout" \
 "$HIST/frozen/F3B3A_IMPLEMENTATION_LOCK.json" \
 "$HERE/work/historical_replay" 32 > "$HERE/logs/historical_replay.log" 2>&1
python - "$HIST" "$HERE/work/historical_replay" <<'PY'
import hashlib,sys
from pathlib import Path
H,R=map(Path,sys.argv[1:3])
def h(p):
 x=hashlib.sha256();x.update(open(p,'rb').read());return x.hexdigest()
for n in ['F3B3A_CANDIDATE_CHARACTERIZATION.tsv','F3B3A_FAMILY_PRIORITY.tsv','F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv','F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.faa','STAGE05F3B3A_GATE.json']:
 a=H/'outputs'/n;b=R/n
 if h(a)!=h(b):raise SystemExit('HOLD historical replay differs '+n)
print('PASS_STAGE04B9_HISTORICAL_F3B3A_EXACT_REPLAY')
PY
python "$HERE/scripts/historical_semantics_audit.py" \
 --r4-map "$R4/outputs/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv" \
 --r3-domtblout "$R3/raw/F3B1R3_RVT_ALL.domtblout" \
 --r3-pair-hits "$R3/outputs/STAGE05F3B1R3_INSTANCE_RVT_PAIR_HITS.tsv.gz" \
 --out "$HERE/outputs/HISTORICAL_MODEL_SEMANTICS_AUDIT.json" | tee "$HERE/logs/historical_semantics.log"
python "$HERE/scripts/strict_mirror_prepare.py" --project-root "$P" --home "$HERE" --semantics-audit "$HERE/outputs/HISTORICAL_MODEL_SEMANTICS_AUDIT.json" | tee "$HERE/logs/strict_mirror_prepare.log"
python "$HERE/scripts/prepare_arch5.py" --project-root "$P" --home "$HERE" | tee "$HERE/logs/arch5_prepare.log"
echo PASS_STAGE04B9_PREPARE_BOTH_BRANCHES_COMPLETE
