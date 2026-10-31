#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_PROJECT_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
B=${F3B4E_BASE_DIR:-$P/ART_MATH_DISCOVERY_STAGE05F3B4E_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_v1_0_0_20261011}
R3=${F3B4E_R3_DIR:-$P/ART_MATH_DISCOVERY_STAGE05F3B4E_G1_RID_RESUME_HOTFIX_v1_0_3_20261012}
mkdir -p "$S"/{outputs,logs,tmp,frozen}
for x in python sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for f in "$B/inputs/F3B4E_RT5.faa" "$B/inputs/F3B4E_PARTNER5.faa" "$B/inputs/F3B4E_ARCHITECTURE5.fna" "$B/outputs/F3B4E_TARGETS_5.tsv" "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json"; do [[ -s "$f" ]] || { echo "[FATAL] missing $f" >&2; exit 2; }; done
sha_eq(){ local f=$1 exp=$2 got; got=$(sha256sum "$f"|awk '{print $1}'); [[ "$got" == "$exp" ]] || { echo "[FATAL] SHA mismatch $f" >&2; exit 2; }; }
sha_eq "$B/inputs/F3B4E_RT5.faa" 2bc4b561e4df459c5ca6fd380054076b3f0a1e5330305f32072d79df565c81bb
sha_eq "$B/inputs/F3B4E_PARTNER5.faa" ea904cfd6fb2bc54bcb7aee1a23ecb379d3a563cce486670e01f807033ff15f6
sha_eq "$B/inputs/F3B4E_ARCHITECTURE5.fna" c496dca26651abd3f83b74ece325cb947d315a5fc11c2c880f932a1cfad75cbe
sha_eq "$B/outputs/F3B4E_TARGETS_5.tsv" 9a631731f1452978ae0c5731e6393f880fae0e820ec12aa0c04c6da1c6530b27
sha_eq "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" 812b34b6d103da1697f5a0e236c24a225604a8adff3192d971de4bba255a2a00
python - "$S/frozen/F3B4E_G1_R4_SINGLEFILE_XML2_LOCK.json" "$S/scripts/ncbi_blast_urlapi_singlefile.py" "$S/scripts/summarize_global.py" <<'PY'
import hashlib,json,sys
from pathlib import Path
out,client,summarizer=sys.argv[1:]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
z={
 'stage':'Stage05F3B4E-G1-R4-single-file-XML2-retrieval-freeze',
 'gate':'PASS_F3B4E_G1_R4_SINGLE_FILE_XML2_RETRIEVAL_HOTFIX_FROZEN_PRE_RESULTS',
 'gate_pass':True,
 'hotfix_role':'RESULT_BLIND_RESULT_RETRIEVAL_FORMAT_ONLY',
 'trigger':'Existing RT5 RID became READY, but R3 repeatedly requested FORMAT_TYPE=XML2 (multiple-file XML2) and then searched the returned payload for single-file XML tags, causing a READY loop.',
 'existing_ready_rid_reuse_first':True,
 'retrieval_format_old':'XML2 (multiple-file report)',
 'retrieval_format_new':'XML2_S (single-file report)',
 'query_membership_changed':False,'query_sequences_changed':False,'database_scope_changed':False,
 'evalue_changed':False,'hitlist_size_changed':False,'filtering_changed':False,'novelty_thresholds_changed':False,
 'resubmit_existing_ready_rid':False,
 'ready_retrieval_retry_n':5,'ready_retrieval_retry_seconds':30,
 'urlapi_client_sha256':sha(client),'summarizer_sha256':sha(summarizer),
 'claim_boundary':'Retrieval-format repair only. READY indicates completed server search; changing XML2 to XML2_S changes only serialization/retrieval, not BLAST search or biological results.'}
Path(out).write_text(json.dumps(z,indent=2)+'\n');print(json.dumps(z,indent=2))
PY

RT_RESUME_ARGS=()
if [[ -n "${F3B4E_RT5_RID:-}" ]]; then RT_RESUME_ARGS=(--resume-rid "$F3B4E_RT5_RID");
elif [[ -s "$R3/outputs/F3B4E_RT5_NCBI_NR.rid.tsv" ]]; then RT_RESUME_ARGS=(--resume-rid-file "$R3/outputs/F3B4E_RT5_NCBI_NR.rid.tsv");
else echo '[FATAL] RT5 READY RID not supplied/found. Set F3B4E_RT5_RID.' >&2; exit 2; fi

echo '[1/3] RT5 vs NCBI nr: reuse existing RID; retrieve SINGLE-FILE XML2'
python "$S/scripts/ncbi_blast_urlapi_singlefile.py" --program blastp --database nr --query "$B/inputs/F3B4E_RT5.faa" --expect 1e-3 --hitlist 200 --filter T --out "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" --xml-out "$S/outputs/F3B4E_RT5_NCBI_NR.xml" --rid-out "$S/outputs/F3B4E_RT5_NCBI_NR.rid.tsv" --max-wait 64800 "${RT_RESUME_ARGS[@]}" 2>&1 | tee "$S/logs/rt5_urlapi_singlefile.log"

echo '[2/3] partner5 vs NCBI nr: submit/resume; retrieve SINGLE-FILE XML2'
PART_ARGS=(); [[ -s "$S/outputs/F3B4E_PARTNER5_NCBI_NR.rid.tsv" ]] && PART_ARGS=(--resume-rid-file "$S/outputs/F3B4E_PARTNER5_NCBI_NR.rid.tsv")
python "$S/scripts/ncbi_blast_urlapi_singlefile.py" --program blastp --database nr --query "$B/inputs/F3B4E_PARTNER5.faa" --expect 1e-3 --hitlist 200 --filter T --out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" --xml-out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.xml" --rid-out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.rid.tsv" --max-wait 64800 "${PART_ARGS[@]}" 2>&1 | tee "$S/logs/partner5_urlapi_singlefile.log"

echo '[3/3] architecture5 vs NCBI nt: submit/resume; retrieve SINGLE-FILE XML2'
ARCH_ARGS=(); [[ -s "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.rid.tsv" ]] && ARCH_ARGS=(--resume-rid-file "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.rid.tsv")
python "$S/scripts/ncbi_blast_urlapi_singlefile.py" --program blastn --database nt --megablast --no-db-override --query "$B/inputs/F3B4E_ARCHITECTURE5.fna" --expect 1e-10 --hitlist 200 --filter T --out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" --xml-out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.xml" --rid-out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.rid.tsv" --max-wait 64800 "${ARCH_ARGS[@]}" 2>&1 | tee "$S/logs/architecture5_urlapi_singlefile.log"
python - "$S/outputs/F3B4E_GLOBAL_URLAPI_RUNTIME.json" <<'PY'
import json,datetime,sys
p=sys.argv[1];json.dump({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'transport':'NCBI Common URL API HTTPS with persistent RID and XML2_S single-file retrieval','database_scope':['nr','nt']},open(p,'w'),indent=2);open(p,'a').write('\n')
PY
PYTHONPATH="$S/scripts" python "$S/scripts/summarize_global.py" "$B/outputs/F3B4E_TARGETS_5.tsv" "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" "$S/outputs/F3B4E_GLOBAL_URLAPI_RUNTIME.json" | tee "$S/logs/summarize_global_urlapi_singlefile.log"
sha256sum "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" "$S/outputs/STAGE05F3B4E_G1_GATE.json" "$S/frozen/F3B4E_G1_R4_SINGLEFILE_XML2_LOCK.json" > "$S/outputs/F3B4E_G1_R4_SHA256.tsv"
echo PASS_STAGE05F3B4E_G1_R4_SINGLE_FILE_XML2_COMPLETE
