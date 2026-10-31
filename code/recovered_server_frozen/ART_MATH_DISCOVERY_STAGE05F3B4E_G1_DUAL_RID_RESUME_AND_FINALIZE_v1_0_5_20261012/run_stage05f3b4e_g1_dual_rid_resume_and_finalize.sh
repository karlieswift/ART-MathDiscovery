#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_PROJECT_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
B=${F3B4E_BASE_DIR:-$P/ART_MATH_DISCOVERY_STAGE05F3B4E_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_v1_0_0_20261011}
RT_RID=${F3B4E_RT5_RID:-CR73H1ST014}
PART_RID=${F3B4E_PARTNER5_RID:-CSE25BRU016}
ARCH_RID=${F3B4E_ARCHITECTURE5_RID:-CSRWCA1M016}
MAX_WAIT=${F3B4E_RESUME_MAX_WAIT_SECONDS:-86400}
mkdir -p "$S"/{outputs,logs,tmp,frozen}
for x in python sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for f in \
  "$B/inputs/F3B4E_RT5.faa" \
  "$B/inputs/F3B4E_PARTNER5.faa" \
  "$B/inputs/F3B4E_ARCHITECTURE5.fna" \
  "$B/outputs/F3B4E_TARGETS_5.tsv" \
  "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" \
  "$B/scripts/finalize.py" \
  "$B/finalize_stage05f3b4e.sh" \
  "$B/collect_stage05f3b4e.sh"
do
  [[ -s "$f" ]] || { echo "[FATAL] missing $f" >&2; exit 2; }
done
sha_eq(){ local f=$1 exp=$2 got; got=$(sha256sum "$f"|awk '{print $1}'); [[ "$got" == "$exp" ]] || { echo "[FATAL] SHA mismatch $f" >&2; exit 2; }; }
sha_eq "$B/inputs/F3B4E_RT5.faa" 2bc4b561e4df459c5ca6fd380054076b3f0a1e5330305f32072d79df565c81bb
sha_eq "$B/inputs/F3B4E_PARTNER5.faa" ea904cfd6fb2bc54bcb7aee1a23ecb379d3a563cce486670e01f807033ff15f6
sha_eq "$B/inputs/F3B4E_ARCHITECTURE5.fna" c496dca26651abd3f83b74ece325cb947d315a5fc11c2c880f932a1cfad75cbe
sha_eq "$B/outputs/F3B4E_TARGETS_5.tsv" 9a631731f1452978ae0c5731e6393f880fae0e820ec12aa0c04c6da1c6530b27
sha_eq "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" 812b34b6d103da1697f5a0e236c24a225604a8adff3192d971de4bba255a2a00

python - "$S/frozen/F3B4E_G1_R5_DUAL_RID_RESUME_LOCK.json" "$S/scripts/ncbi_blast_urlapi_singlefile.py" "$S/scripts/summarize_global.py" "$RT_RID" "$PART_RID" "$ARCH_RID" <<'PY'
import hashlib,json,sys
from pathlib import Path
out,client,summarizer,rt,partner,arch=sys.argv[1:]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
z={
 'stage':'Stage05F3B4E-G1-R5-dual-RID-resume-and-finalize',
 'gate':'PASS_F3B4E_G1_R5_EXISTING_RIDS_REUSED_WITHOUT_RESUBMISSION',
 'gate_pass':True,
 'hotfix_role':'OPERATIONAL_CONSOLIDATION_OF_ALREADY_SUBMITTED_NCBI_RIDS_ONLY',
 'rt5_rid':rt,'partner5_rid':partner,'architecture5_rid':arch,
 'submission_policy':'NO_NEW_SUBMISSION; existing RIDs only',
 'query_membership_changed':False,'query_sequences_changed':False,'database_scope_changed':False,
 'evalue_changed':False,'hitlist_size_changed':False,'filtering_changed':False,'novelty_thresholds_changed':False,
 'rt5_search':'blastp vs nr; E<=1e-3; hitlist 200; filter T',
 'partner5_search':'blastp vs nr; E<=1e-3; hitlist 200; filter T',
 'architecture5_search':'megablast vs nt; E<=1e-10; hitlist 200; filter T; NO_DATABASE_OVERRIDE',
 'urlapi_client_sha256':sha(client),'summarizer_sha256':sha(summarizer),
 'claim_boundary':'Operational resume/finalization only. The three NCBI searches were already submitted under frozen F3B4E rules; this runner reuses their RIDs, retrieves results, summarizes, copies G1 outputs into the canonical base stage, and runs the pre-existing finalizer.'}
Path(out).write_text(json.dumps(z,indent=2)+'\n')
print(json.dumps(z,indent=2))
PY

echo "=== EXISTING NCBI RID STATE ==="
echo "RT5          $RT_RID"
echo "PARTNER5     $PART_RID"
echo "ARCHITECTURE $ARCH_RID"
echo "No RID will be resubmitted by this runner."

echo "[1/4] retrieve/resume RT5 existing RID"
python "$S/scripts/ncbi_blast_urlapi_singlefile.py" \
  --program blastp --database nr \
  --query "$B/inputs/F3B4E_RT5.faa" \
  --expect 1e-3 --hitlist 200 --filter T \
  --out "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" \
  --xml-out "$S/outputs/F3B4E_RT5_NCBI_NR.xml" \
  --rid-out "$S/outputs/F3B4E_RT5_NCBI_NR.rid.tsv" \
  --resume-rid "$RT_RID" --max-wait "$MAX_WAIT" \
  > "$S/logs/rt5_resume.log" 2>&1
 tail -5 "$S/logs/rt5_resume.log"

echo "[2/4] resume PARTNER5 and ARCHITECTURE5 in parallel; existing RIDs only"
set +e
python "$S/scripts/ncbi_blast_urlapi_singlefile.py" \
  --program blastp --database nr \
  --query "$B/inputs/F3B4E_PARTNER5.faa" \
  --expect 1e-3 --hitlist 200 --filter T \
  --out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" \
  --xml-out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.xml" \
  --rid-out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.rid.tsv" \
  --resume-rid "$PART_RID" --max-wait "$MAX_WAIT" \
  > "$S/logs/partner5_resume.log" 2>&1 &
PP=$!
python "$S/scripts/ncbi_blast_urlapi_singlefile.py" \
  --program blastn --database nt --megablast --no-db-override \
  --query "$B/inputs/F3B4E_ARCHITECTURE5.fna" \
  --expect 1e-10 --hitlist 200 --filter T \
  --out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" \
  --xml-out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.xml" \
  --rid-out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.rid.tsv" \
  --resume-rid "$ARCH_RID" --max-wait "$MAX_WAIT" \
  > "$S/logs/architecture5_resume.log" 2>&1 &
AP=$!
set -e
while kill -0 "$PP" 2>/dev/null || kill -0 "$AP" 2>/dev/null; do
  echo "--- $(date '+%Y-%m-%d %H:%M:%S') ---"
  echo "PARTNER5:     $(tail -1 "$S/logs/partner5_resume.log" 2>/dev/null || true)"
  echo "ARCHITECTURE5:$(tail -1 "$S/logs/architecture5_resume.log" 2>/dev/null || true)"
  sleep 600
done
set +e
wait "$PP"; PRC=$?
wait "$AP"; ARC=$?
set -e
if [[ "$PRC" -ne 0 || "$ARC" -ne 0 ]]; then
  echo "[FATAL] resume incomplete: partner_rc=$PRC architecture_rc=$ARC" >&2
  echo "=== partner tail ===" >&2; tail -20 "$S/logs/partner5_resume.log" >&2 || true
  echo "=== architecture tail ===" >&2; tail -20 "$S/logs/architecture5_resume.log" >&2 || true
  echo "RIDs are preserved. Re-run this same command later; do not submit new searches." >&2
  exit 3
fi

echo "[3/4] summarize frozen global public audit"
python - "$S/outputs/F3B4E_GLOBAL_URLAPI_RUNTIME.json" "$RT_RID" "$PART_RID" "$ARCH_RID" <<'PY'
import json,datetime,sys
p,rt,pa,ar=sys.argv[1:]
z={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
   'transport':'NCBI Common URL API HTTPS; existing RID resume; XML2_S retrieval',
   'database_scope':['nr','nt'],'rt5_rid':rt,'partner5_rid':pa,'architecture5_rid':ar,
   'new_submission_performed':False}
open(p,'w').write(json.dumps(z,indent=2)+'\n')
PY
PYTHONPATH="$S/scripts" python "$S/scripts/summarize_global.py" \
  "$B/outputs/F3B4E_TARGETS_5.tsv" \
  "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" \
  "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" \
  "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" \
  "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" \
  "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" \
  "$S/outputs/F3B4E_GLOBAL_URLAPI_RUNTIME.json" \
  | tee "$S/logs/summarize_global.log"

for f in \
  F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv \
  STAGE05F3B4E_G1_GATE.json \
  F3B4E_RT5_NCBI_NR.tsv \
  F3B4E_PARTNER5_NCBI_NR.tsv \
  F3B4E_ARCHITECTURE5_NCBI_NT.tsv \
  F3B4E_GLOBAL_URLAPI_RUNTIME.json
do
  cp -f "$S/outputs/$f" "$B/outputs/$f"
done
cp -f "$S/frozen/F3B4E_G1_R5_DUAL_RID_RESUME_LOCK.json" "$B/outputs/F3B4E_G1_R5_DUAL_RID_RESUME_LOCK.json"

[[ -s "$B/outputs/STAGE05F3B4E_G2_GATE.json" ]] || { echo '[FATAL] Public ART8 G2 gate missing in canonical base. Run the already-completed Public ART8 R1 on this shared project root or copy its canonical outputs.' >&2; exit 4; }
[[ -s "$B/outputs/F3B4E_PUBLIC_ART8_FAMILY_SUMMARY_5.tsv" ]] || { echo '[FATAL] Public ART8 family summary missing in canonical base.' >&2; exit 4; }

echo "[4/4] canonical F3B4E finalization + collection"
cd "$B"
./finalize_stage05f3b4e.sh
./collect_stage05f3b4e.sh

echo
cat outputs/STAGE05F3B4E_GATE.json

echo
ls -lh "$P"/ART_STAGE05F3B4E_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_RESULTS_*.tar.gz 2>/dev/null | tail -5 || true

echo PASS_STAGE05F3B4E_G1_R5_DUAL_RID_RESUME_AND_FULL_FINALIZATION_COMPLETE
