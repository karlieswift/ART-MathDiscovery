#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_PROJECT_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
B=${F3B4E_BASE_DIR:-$P/ART_MATH_DISCOVERY_STAGE05F3B4E_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_v1_0_0_20261011}
R2=${F3B4E_R2_DIR:-$P/ART_MATH_DISCOVERY_STAGE05F3B4E_G1_HTTP_URLAPI_HOTFIX_v1_0_2_20261011}
mkdir -p "$S"/{outputs,logs,tmp,frozen}
for x in python sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for f in "$B/inputs/F3B4E_RT5.faa" "$B/inputs/F3B4E_PARTNER5.faa" "$B/inputs/F3B4E_ARCHITECTURE5.fna" "$B/outputs/F3B4E_TARGETS_5.tsv" "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json"; do [[ -s "$f" ]] || { echo "[FATAL] missing $f" >&2; exit 2; }; done
sha_eq(){ local f=$1 exp=$2 got; got=$(sha256sum "$f"|awk '{print $1}'); [[ "$got" == "$exp" ]] || { echo "[FATAL] SHA mismatch $f" >&2; exit 2; }; }
sha_eq "$B/inputs/F3B4E_RT5.faa" 2bc4b561e4df459c5ca6fd380054076b3f0a1e5330305f32072d79df565c81bb
sha_eq "$B/inputs/F3B4E_PARTNER5.faa" ea904cfd6fb2bc54bcb7aee1a23ecb379d3a563cce486670e01f807033ff15f6
sha_eq "$B/inputs/F3B4E_ARCHITECTURE5.fna" c496dca26651abd3f83b74ece325cb947d315a5fc11c2c880f932a1cfad75cbe
sha_eq "$B/outputs/F3B4E_TARGETS_5.tsv" 9a631731f1452978ae0c5731e6393f880fae0e820ec12aa0c04c6da1c6530b27
sha_eq "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" 812b34b6d103da1697f5a0e236c24a225604a8adff3192d971de4bba255a2a00
python - "$S/frozen/F3B4E_G1_R3_RID_RESUME_LOCK.json" "$S/scripts/ncbi_blast_urlapi_resume.py" "$S/scripts/summarize_global.py" <<'PY'
import hashlib,json,sys
from pathlib import Path
out,client,summarizer=sys.argv[1:]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
z={'stage':'Stage05F3B4E-G1-R3-RID-resume-freeze','gate':'PASS_F3B4E_G1_R3_EXISTING_RID_RESUME_POLICY_FROZEN_PRE_RESULTS','gate_pass':True,'hotfix_role':'RESULT_BLIND_TRANSPORT_STATE_RESUME_ONLY','trigger':'R2 obtained RT5 RID CR73H1ST014 but fixed 10800-second client wait expired while server still reported pending','query_membership_changed':False,'query_sequences_changed':False,'database_scope_changed':False,'evalue_changed':False,'hitlist_size_changed':False,'filtering_changed':False,'novelty_thresholds_changed':False,'existing_rid_reuse_first':True,'poll_interval_seconds':60,'resume_invocation_max_wait_seconds':64800,'timeout_semantics':'client timeout only; preserve RID and resume again; never interpret timeout as no-hit','resubmit_policy':'do not resubmit a valid WAITING/READY RID; resubmit only after explicit expired/UNKNOWN/FAILED state and no retrievable result','urlapi_client_sha256':sha(client),'summarizer_sha256':sha(summarizer),'claim_boundary':'Transport-state persistence only; scientific searches remain RT5 blastp/nr, partner5 blastp/nr, architecture5 megablast/nt under the frozen F3B4E rules.'}
Path(out).write_text(json.dumps(z,indent=2)+'\n');print(json.dumps(z,indent=2))
PY

RT_OLD="$R2/outputs/F3B4E_RT5_NCBI_NR.rid.tsv"
RT_RESUME_ARGS=()
if [[ -s "$RT_OLD" ]]; then RT_RESUME_ARGS=(--resume-rid-file "$RT_OLD"); fi
if [[ -n "${F3B4E_RT5_RID:-}" ]]; then RT_RESUME_ARGS=(--resume-rid "$F3B4E_RT5_RID"); fi

echo '[1/3] RT5 vs NCBI nr: resume existing RID first'
python "$S/scripts/ncbi_blast_urlapi_resume.py" --program blastp --database nr --query "$B/inputs/F3B4E_RT5.faa" --expect 1e-3 --hitlist 200 --filter T --out "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" --xml-out "$S/outputs/F3B4E_RT5_NCBI_NR.xml" --rid-out "$S/outputs/F3B4E_RT5_NCBI_NR.rid.tsv" --max-wait 64800 "${RT_RESUME_ARGS[@]}" 2>&1 | tee "$S/logs/rt5_urlapi_resume.log"

echo '[2/3] partner5 vs NCBI nr: new or resume local R3 RID'
PART_ARGS=(); [[ -s "$S/outputs/F3B4E_PARTNER5_NCBI_NR.rid.tsv" ]] && PART_ARGS=(--resume-rid-file "$S/outputs/F3B4E_PARTNER5_NCBI_NR.rid.tsv")
python "$S/scripts/ncbi_blast_urlapi_resume.py" --program blastp --database nr --query "$B/inputs/F3B4E_PARTNER5.faa" --expect 1e-3 --hitlist 200 --filter T --out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" --xml-out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.xml" --rid-out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.rid.tsv" --max-wait 64800 "${PART_ARGS[@]}" 2>&1 | tee "$S/logs/partner5_urlapi_resume.log"

echo '[3/3] architecture5 vs NCBI nt: new or resume local R3 RID'
ARCH_ARGS=(); [[ -s "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.rid.tsv" ]] && ARCH_ARGS=(--resume-rid-file "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.rid.tsv")
python "$S/scripts/ncbi_blast_urlapi_resume.py" --program blastn --database nt --megablast --no-db-override --query "$B/inputs/F3B4E_ARCHITECTURE5.fna" --expect 1e-10 --hitlist 200 --filter T --out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" --xml-out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.xml" --rid-out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.rid.tsv" --max-wait 64800 "${ARCH_ARGS[@]}" 2>&1 | tee "$S/logs/architecture5_urlapi_resume.log"
python - "$S/outputs/F3B4E_GLOBAL_URLAPI_RUNTIME.json" <<'PY'
import json,datetime,sys
p=sys.argv[1];json.dump({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'transport':'NCBI Common URL API HTTPS with persistent RID resume','database_scope':['nr','nt']},open(p,'w'),indent=2);open(p,'a').write('\n')
PY
PYTHONPATH="$S/scripts" python "$S/scripts/summarize_global.py" "$B/outputs/F3B4E_TARGETS_5.tsv" "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" "$S/outputs/F3B4E_GLOBAL_URLAPI_RUNTIME.json" | tee "$S/logs/summarize_global_urlapi_resume.log"
sha256sum "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" "$S/outputs/STAGE05F3B4E_G1_GATE.json" "$S/frozen/F3B4E_G1_R3_RID_RESUME_LOCK.json" > "$S/outputs/F3B4E_G1_R3_SHA256.tsv"
echo PASS_STAGE05F3B4E_G1_R3_RID_RESUME_COMPLETE
