#!/usr/bin/env bash
set -euo pipefail
P=${1:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
HERE=$(cd "$(dirname "$0")" && pwd)
HIST="$P/ART_MATH_DISCOVERY_STAGE05F3B3A_POSTFREEZE_PRIORITY_AND_CATALYTIC_TOPOLOGY_v1_0_0_20261011"
R3="$P/ART_MATH_DISCOVERY_STAGE05F3B1R3_QUERY_NAMESPACE_RESCUE_AND_F3B2_FAMILY_COMPRESSION_v1_0_0_20261011"
R4="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_v1_0_0_20261011"
D="$P/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_DISTRIBUTED_ACQUISITION_R2_v1_0_2_20261013"
A="$P/ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04A_60G_SATURATION_CENSUS_v1_0_0_20261016"
B="$P/ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B1_DELTA_FAMILY_REPLICATION_CENSUS_v1_0_0_20261016"
OUT="$HERE/outputs"; WORK="$HERE/work"; FROZEN="$HERE/frozen"; LOG="$HERE/logs"
mkdir -p "$OUT" "$WORK/historical_replay" "$FROZEN" "$LOG"

"$HERE/tests/run_selftest.sh"

HISTSCRIPT="$HIST/scripts/stage05f3b3a.py"
EXACT="$HERE/frozen/historical_stage05f3b3a_exact.py"
EXPECTED_SCRIPT_SHA=d222f638c627566909765012874c7aaa531a1a79b5e67a0d73f726f184a30ed8
[[ -s "$HISTSCRIPT" && -s "$HIST/frozen/F3B3A_IMPLEMENTATION_LOCK.json" && -s "$HIST/outputs/F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv" ]] || { echo '[HOLD] historical F3B3A runtime/output missing' >&2; exit 2; }
[[ $(sha256sum "$HISTSCRIPT"|awk '{print $1}') == "$EXPECTED_SCRIPT_SHA" ]] || { echo '[HOLD] historical F3B3A script SHA mismatch' >&2; exit 2; }
[[ $(sha256sum "$EXACT"|awk '{print $1}') == "$EXPECTED_SCRIPT_SHA" ]] || { echo '[HOLD] bundled historical F3B3A exact script SHA mismatch' >&2; exit 2; }

# Freeze rule + input hashes before delta characterization.
python - "$HERE" "$HIST" "$R3" "$R4" "$D" "$A" "$B" "$FROZEN/STAGE04B2_RULE_AND_INPUT_LOCK.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
HERE,HIST,R3,R4,D,A,B,out=map(Path,sys.argv[1:9])
def h(p):
 x=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):x.update(b)
 return x.hexdigest()
paths={
'historical_script':HIST/'scripts/stage05f3b3a.py',
'historical_lock':HIST/'frozen/F3B3A_IMPLEMENTATION_LOCK.json',
'historical_shortlist':HIST/'outputs/F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv',
'r4_map':R4/'outputs/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv',
'r3_pair':R3/'outputs/STAGE05F3B1R3_INSTANCE_RVT_PAIR_HITS.tsv.gz',
'r3_dom':R3/'raw/F3B1R3_RVT_ALL.domtblout',
'delta_gate':D/'outputs/STAGE05F3V2_GATE.json',
'delta_formal_instances':D/'outputs/F3V2_FORMAL_INSTANCES.tsv',
'delta_manifest':D/'frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv',
'stage04a_gate':A/'outputs/STAGE04A_60G_SATURATION_CENSUS.json',
'stage04a_membership':A/'outputs/COMBINED_FAMILY_MEMBERSHIP.tsv',
'stage04b1_gate':B/'outputs/STAGE04B1_DELTA_FAMILY_REPLICATION_CENSUS.json',
'projection_script':HERE/'scripts/delta_priority_projection.py',
'semantics_script':HERE/'scripts/historical_semantics_audit.py',
'rule_core':HERE/'scripts/rule_core.py',
'bundled_historical_script':HERE/'frozen/historical_stage05f3b3a_exact.py',
}
for k,p in paths.items():
 if not p.is_file(): raise SystemExit(f'HOLD missing input {k}: {p}')
g={
 'stage':'Stage05F3V2-Stage04B2-rule-input-freeze',
 'gate':'PASS_F3V2_STAGE04B2_HISTORICAL_RULE_AND_DELTA_INPUTS_FROZEN_PRE_PROJECTION',
 'gate_pass':True,
 'hashes':{k:h(p) for k,p in paths.items()},
 'historical_rule_source_date':'2026-10-11',
 'historical_rule':'four axes x8: replication, max ORF length, formal-model coverage, motif-C topology; preserve axis order; deterministic backfill to 32',
 'motif_c_like_regex':'[YF]xDD',
 'family_universe':'Stage04A DELTA_ONLY combined families only',
 'known_rt_annotation_used_for_ranking':False,'novelty_used_for_ranking':False,'structure_used_for_ranking':False,
 'retuning_authorized':False,'science_change':False,
 'claim_boundary':'The historical F3B3A priority implementation predates the 54-Gb delta. This lock fixes the exact source implementation and delta-only combined-family input universe before Stage04B2 projection results are generated.'
}
with open(out,'w') as f:json.dump(g,f,indent=2);f.write('\n')
print(json.dumps(g,indent=2))
PY

# Exact historical replay into this package; never modify historical outputs.
rm -rf "$WORK/historical_replay"; mkdir -p "$WORK/historical_replay"
python "$EXACT" \
  "$R4/outputs/STAGE05F3B2_R4_GATE.json" \
  "$R4/outputs/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv" \
  "$R4/outputs/F3B2_R4_FAMILY_SUMMARY_WITH_SOURCE.tsv" \
  "$R3/outputs/STAGE05F3B1R3_FORMAL_INSTANCE_IDS.txt" \
  "$R3/raw/F3B1R3_NAMESPACE_INSTANCES.faa" \
  "$R3/outputs/STAGE05F3B1R3_INSTANCE_RVT_PAIR_HITS.tsv.gz" \
  "$R3/raw/F3B1R3_RVT_ALL.domtblout" \
  "$HIST/frozen/F3B3A_IMPLEMENTATION_LOCK.json" \
  "$WORK/historical_replay" 32 \
  > "$LOG/historical_exact_replay.log" 2>&1

python - "$HIST" "$WORK/historical_replay" "$OUT/HISTORICAL_F3B3A_REPLAY_GATE.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
H,R,O=map(Path,sys.argv[1:4])
def h(p):
 x=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):x.update(b)
 return x.hexdigest()
files=['F3B3A_CANDIDATE_CHARACTERIZATION.tsv','F3B3A_FAMILY_PRIORITY.tsv','F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv','F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.faa','STAGE05F3B3A_GATE.json']
rows=[];ok=True
for n in files:
 a=H/'outputs'/n;b=R/n
 if not a.is_file() or not b.is_file(): raise SystemExit('HOLD missing historical replay comparison '+n)
 x=h(a);y=h(b);same=x==y;ok &= same;rows.append({'file':n,'historical_sha256':x,'replay_sha256':y,'byte_identical':same})
g={'stage':'Stage05F3V2-Stage04B2-historical-replay','gate':'PASS_HISTORICAL_F3B3A_EXACT_OUTPUT_REPLAY' if ok else 'HOLD_HISTORICAL_F3B3A_REPLAY_DIFFERENCE','gate_pass':bool(ok),'comparisons':rows,'retuning_authorized':False,'science_change':False}
with open(O,'w') as f:json.dump(g,f,indent=2);f.write('\n')
print(json.dumps(g,indent=2))
if not ok: raise SystemExit(2)
PY

# Recover the exact model-count meaning from historical R3/R4 evidence.
python "$HERE/scripts/historical_semantics_audit.py" \
  --r4-map "$R4/outputs/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv" \
  --r3-domtblout "$R3/raw/F3B1R3_RVT_ALL.domtblout" \
  --r3-pair-hits "$R3/outputs/STAGE05F3B1R3_INSTANCE_RVT_PAIR_HITS.tsv.gz" \
  --out "$OUT/HISTORICAL_MODEL_SEMANTICS_AUDIT.json" \
  | tee "$LOG/historical_model_semantics.log"

# Apply the exact historical ranking logic to the frozen DELTA_ONLY combined-family universe.
python "$HERE/scripts/delta_priority_projection.py" \
  --stage04a-root "$A" \
  --stage04b1-root "$B" \
  --delta-root "$D" \
  --semantics-audit "$OUT/HISTORICAL_MODEL_SEMANTICS_AUDIT.json" \
  --historical-replay-gate "$OUT/HISTORICAL_F3B3A_REPLAY_GATE.json" \
  --out "$OUT" \
  | tee "$LOG/delta_exact_priority_projection.log"

sha256sum \
  "$FROZEN/STAGE04B2_RULE_AND_INPUT_LOCK.json" \
  "$OUT/HISTORICAL_F3B3A_REPLAY_GATE.json" \
  "$OUT/HISTORICAL_MODEL_SEMANTICS_AUDIT.json" \
  "$OUT/STAGE04B2_GATE.json" \
  "$OUT/F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.tsv" \
  "$OUT/F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.faa" \
  > "$OUT/STAGE04B2_PRINCIPAL_SHA256.txt"

echo
echo '=== FINAL GATE ==='
cat "$OUT/STAGE04B2_GATE.json"
echo
echo '=== SHORTLIST 32 ==='
column -t -s $'\t' "$OUT/F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.tsv" | head -33 || true
echo
echo PASS_STAGE05F3V2_STAGE04B2_HISTORICAL_F3B3A_EXACT_PRIORITY_PROJECTION_COMPLETE
