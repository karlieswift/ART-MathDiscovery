#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_PROJECT_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
B=${F3B4E_BASE_DIR:-$P/ART_MATH_DISCOVERY_STAGE05F3B4E_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_v1_0_0_20261011}
mkdir -p "$S"/{outputs,logs,tmp,frozen}
for x in python sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for f in "$B/inputs/F3B4E_RT5.faa" "$B/inputs/F3B4E_PARTNER5.faa" "$B/inputs/F3B4E_ARCHITECTURE5.fna" "$B/outputs/F3B4E_TARGETS_5.tsv" "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json"; do [[ -s "$f" ]] || { echo "[FATAL] missing prepared base input: $f" >&2; exit 2; }; done
sha_eq(){ local f=$1 exp=$2; local got; got=$(sha256sum "$f" | awk '{print $1}'); [[ "$got" == "$exp" ]] || { echo "[FATAL] SHA mismatch $f got=$got expected=$exp" >&2; exit 2; }; }
sha_eq "$B/inputs/F3B4E_RT5.faa" 2bc4b561e4df459c5ca6fd380054076b3f0a1e5330305f32072d79df565c81bb
sha_eq "$B/inputs/F3B4E_PARTNER5.faa" ea904cfd6fb2bc54bcb7aee1a23ecb379d3a563cce486670e01f807033ff15f6
sha_eq "$B/inputs/F3B4E_ARCHITECTURE5.fna" c496dca26651abd3f83b74ece325cb947d315a5fc11c2c880f932a1cfad75cbe
sha_eq "$B/outputs/F3B4E_TARGETS_5.tsv" 9a631731f1452978ae0c5731e6393f880fae0e820ec12aa0c04c6da1c6530b27
sha_eq "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" 812b34b6d103da1697f5a0e236c24a225604a8adff3192d971de4bba255a2a00
python - "$S/frozen/F3B4E_G1_R2_URLAPI_TRANSPORT_LOCK.json" "$S/scripts/ncbi_blast_urlapi.py" "$S/scripts/summarize_global.py" <<'PY'
import hashlib,json,sys,datetime
from pathlib import Path
out,client,summarizer=sys.argv[1:]
def sha(p):
 h=hashlib.sha256(); h.update(Path(p).read_bytes()); return h.hexdigest()
z={
 'stage':'Stage05F3B4E-G1-R2-transport-freeze','gate':'PASS_F3B4E_G1_R2_HTTPS_COMMON_URL_API_TRANSPORT_FROZEN_PRE_RESULTS','gate_pass':True,
 'hotfix_role':'RESULT_BLIND_TRANSPORT_ONLY',
 'reason':'BLAST+ -remote failed identically on as-4125 and r750-1 with CRPCClient Blast4 request connection stream bad state before any G1 result was produced',
 'old_transport':'NCBI BLAST+ C++ Blast4 RPC via -remote','new_transport':'official NCBI BLAST Common URL API over HTTPS; XML2 result retrieval parsed to the original frozen TSV schema',
 'query_membership_changed':False,'query_sequences_changed':False,'database_scope_changed':False,'evalue_changed':False,'hitlist_size_changed':False,'low_complexity_filtering_changed':False,'novelty_thresholds_changed':False,
 'rt_search':{'program':'blastp','database':'nr','expect':'1e-3','hitlist_size':200,'FILTER':'T'},
 'partner_search':{'program':'blastp','database':'nr','expect':'1e-3','hitlist_size':200,'FILTER':'T'},
 'architecture_search':{'program':'blastn','MEGABLAST':'on','database':'nt','NO_DATABASE_OVERRIDE':'true','expect':'1e-10','hitlist_size':200,'FILTER':'T'},
 'poll_policy':'serial three batched searches; wait >=60 s between polls for a RID; maximum 3 h per search',
 'urlapi_client_sha256':sha(client),'summarizer_sha256':sha(summarizer),
 'claim_boundary':'Transport repair only. HTTPS URL API results are a dated NCBI nr/nt public audit; no-hit is not absolute global novelty.'}
Path(out).write_text(json.dumps(z,indent=2)+'\n')
print(json.dumps(z,indent=2))
PY
echo '[1/3] RT5 vs NCBI nr via HTTPS Common URL API'
python "$S/scripts/ncbi_blast_urlapi.py" --program blastp --database nr --query "$B/inputs/F3B4E_RT5.faa" --expect 1e-3 --hitlist 200 --filter T --out "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" --xml-out "$S/outputs/F3B4E_RT5_NCBI_NR.xml" --rid-out "$S/outputs/F3B4E_RT5_NCBI_NR.rid.tsv" 2>&1 | tee "$S/logs/rt5_urlapi.log"
echo '[2/3] partner5 vs NCBI nr via HTTPS Common URL API'
python "$S/scripts/ncbi_blast_urlapi.py" --program blastp --database nr --query "$B/inputs/F3B4E_PARTNER5.faa" --expect 1e-3 --hitlist 200 --filter T --out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" --xml-out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.xml" --rid-out "$S/outputs/F3B4E_PARTNER5_NCBI_NR.rid.tsv" 2>&1 | tee "$S/logs/partner5_urlapi.log"
echo '[3/3] architecture5 vs NCBI nt via HTTPS Common URL API / MEGABLAST'
python "$S/scripts/ncbi_blast_urlapi.py" --program blastn --database nt --megablast --no-db-override --query "$B/inputs/F3B4E_ARCHITECTURE5.fna" --expect 1e-10 --hitlist 200 --filter T --out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" --xml-out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.xml" --rid-out "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.rid.tsv" 2>&1 | tee "$S/logs/architecture5_urlapi.log"
python - "$S/outputs/F3B4E_GLOBAL_URLAPI_RUNTIME.json" <<'PY'
import json,datetime,sys
json.dump({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'transport':'NCBI Common URL API HTTPS','database_scope':['nr','nt'],'note':'R2 transport hotfix after Blast4 RPC failure'},open(sys.argv[1],'w'),indent=2);open(sys.argv[1],'a').write('\n')
PY
PYTHONPATH="$S/scripts" python "$S/scripts/summarize_global.py" "$B/outputs/F3B4E_TARGETS_5.tsv" "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" "$S/outputs/F3B4E_GLOBAL_URLAPI_RUNTIME.json" | tee "$S/logs/summarize_global_urlapi.log"
sha256sum "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" "$S/outputs/STAGE05F3B4E_G1_GATE.json" "$S/frozen/F3B4E_G1_R2_URLAPI_TRANSPORT_LOCK.json" > "$S/outputs/F3B4E_G1_R2_SHA256.tsv"
echo PASS_STAGE05F3B4E_G1_R2_HTTPS_URLAPI_COMPLETE
