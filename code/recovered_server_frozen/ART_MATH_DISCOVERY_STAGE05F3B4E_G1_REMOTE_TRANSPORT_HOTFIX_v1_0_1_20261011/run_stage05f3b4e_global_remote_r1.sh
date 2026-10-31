#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
BASE="$P/ART_MATH_DISCOVERY_STAGE05F3B4E_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_v1_0_0_20261011"
mkdir -p "$S"/{outputs,logs,work,tmp,frozen}
for x in blastp blastn python sha256sum awk; do
  command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }
done
for f in \
  "$BASE/inputs/F3B4E_RT5.faa" \
  "$BASE/inputs/F3B4E_PARTNER5.faa" \
  "$BASE/inputs/F3B4E_ARCHITECTURE5.fna" \
  "$BASE/outputs/F3B4E_TARGETS_5.tsv" \
  "$BASE/frozen/F3B4E_IMPLEMENTATION_LOCK.json" \
  "$BASE/scripts/summarize_global.py" \
  "$BASE/scripts/common.py"
do
  [[ -s "$f" ]] || { echo "[FATAL] missing base input: $f; run F3B4E prepare first" >&2; exit 2; }
done

# Copy only immutable prepared inputs/scripts needed for transport repair.
cp -f "$BASE/inputs/F3B4E_RT5.faa" "$S/work/RT5.faa"
cp -f "$BASE/inputs/F3B4E_PARTNER5.faa" "$S/work/PARTNER5.faa"
cp -f "$BASE/inputs/F3B4E_ARCHITECTURE5.fna" "$S/work/ARCH5.fna"
cp -f "$BASE/outputs/F3B4E_TARGETS_5.tsv" "$S/work/TARGETS5.tsv"
cp -f "$BASE/frozen/F3B4E_IMPLEMENTATION_LOCK.json" "$S/work/BASE_LOCK.json"
cp -f "$BASE/scripts/summarize_global.py" "$S/scripts/summarize_global.py"
cp -f "$BASE/scripts/common.py" "$S/scripts/common.py"

LOCK="$S/frozen/F3B4E_G1_R1_TRANSPORT_LOCK.json"
if [[ ! -s "$LOCK" ]]; then
python - "$S" "$BASE" "$LOCK" <<'PY'
import hashlib,json,sys,datetime,subprocess
from pathlib import Path
S,BASE,OUT=sys.argv[1:]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def cmd(x):
 try:return subprocess.check_output(x,shell=True,text=True,stderr=subprocess.STDOUT).strip()
 except:return ''
z={
 'stage':'Stage05F3B4E-G1-R1-transport-freeze',
 'gate':'PASS_F3B4E_G1_R1_PER_QUERY_REMOTE_TRANSPORT_HOTFIX_FROZEN_PRE_RESULTS',
 'gate_pass':True,
 'hotfix_role':'RESULT_BLIND_REMOTE_TRANSPORT_ONLY',
 'original_failure':'batched NCBI remote BLASTP failed twice before any RT5 result table was produced',
 'query_membership_changed':False,
 'query_sequences_changed':False,
 'database_scope_changed':False,
 'evalue_changed':False,
 'max_target_seqs_changed':False,
 'seg_dust_changed':False,
 'novelty_thresholds_changed':False,
 'transport_change':'batched five-query request -> one query per request; independent retry/backoff; optional exact same-command local BLASTDB fallback if nr/nt are locally configured',
 'remote_retry_attempts_per_query':6,
 'remote_backoff_seconds':[30,60,120,180,300,480],
 'remote_per_attempt_timeout_seconds':1800,
 'rt_fasta_sha256':sha(Path(S)/'work/RT5.faa'),
 'partner_fasta_sha256':sha(Path(S)/'work/PARTNER5.faa'),
 'architecture_fasta_sha256':sha(Path(S)/'work/ARCH5.fna'),
 'targets_sha256':sha(Path(S)/'work/TARGETS5.tsv'),
 'base_lock_sha256':sha(Path(S)/'work/BASE_LOCK.json'),
 'summarizer_sha256':sha(Path(S)/'scripts/summarize_global.py'),
 'blastp_version':cmd('blastp -version | head -1'),
 'blastn_version':cmd('blastn -version | head -1'),
 'claim_boundary':'Transport repair only. A request failure is technical-censored and is never interpreted as a no-hit or novelty result.'
}
Path(OUT).write_text(json.dumps(z,indent=2)+'\n')
print(json.dumps(z,indent=2))
PY
chmod a-w "$LOCK"
else
  echo '[SKIP] transport lock exists'
  cat "$LOCK"
fi

OUTFMT='6 qseqid saccver pident length qstart qend sstart send evalue bitscore qlen slen qcovs stitle sscinames'
BACKOFF=(30 60 120 180 300 480)

split_fasta(){
  local in=$1 outdir=$2
  rm -rf "$outdir"; mkdir -p "$outdir"
  python - "$in" "$outdir" <<'PY'
import sys,re
from pathlib import Path
inp,out=sys.argv[1:]
name=None; lines=[]; n=0
def flush():
 global name,lines,n
 if name is None:return
 n+=1
 safe=re.sub(r'[^A-Za-z0-9_.-]+','_',name.split()[0])
 Path(out,f'{n:02d}__{safe}.fa').write_text('>'+name+'\n'+''.join(lines)+'\n')
with open(inp) as f:
 for line in f:
  line=line.strip()
  if not line: continue
  if line.startswith('>'):
   flush(); name=line[1:]; lines=[]
  else: lines.append(line)
flush()
print(n)
PY
}

run_one_remote(){
  local kind=$1 query=$2 db=$3 out=$4 log=$5
  local attempt rc
  : > "$log"
  for attempt in 1 2 3 4 5 6; do
    echo "[REMOTE] kind=$kind query=$(basename "$query") db=$db attempt=$attempt" | tee -a "$log"
    set +e
    if [[ "$kind" == protein ]]; then
      timeout 1800 blastp -query "$query" -db "$db" -remote -evalue 1e-3 -max_target_seqs 200 -seg yes -outfmt "$OUTFMT" > "$out.tmp" 2>> "$log"
      rc=$?
    else
      timeout 1800 blastn -task megablast -query "$query" -db "$db" -remote -evalue 1e-10 -max_target_seqs 200 -dust yes -outfmt "$OUTFMT" > "$out.tmp" 2>> "$log"
      rc=$?
    fi
    set -e
    if [[ $rc -eq 0 ]]; then mv -f "$out.tmp" "$out"; return 0; fi
    rm -f "$out.tmp"
    echo "[REMOTE_FAIL] rc=$rc backoff=${BACKOFF[$((attempt-1))]}s" | tee -a "$log"
    sleep "${BACKOFF[$((attempt-1))]}"
  done
  return 1
}

run_one_local_if_available(){
  local kind=$1 query=$2 db=$3 out=$4 log=$5
  command -v blastdbcmd >/dev/null || return 1
  blastdbcmd -db "$db" -info >/dev/null 2>&1 || return 1
  echo "[LOCAL_FALLBACK] query=$(basename "$query") db=$db" | tee -a "$log"
  if [[ "$kind" == protein ]]; then
    blastp -query "$query" -db "$db" -evalue 1e-3 -max_target_seqs 200 -seg yes -outfmt "$OUTFMT" > "$out" 2>> "$log"
  else
    blastn -task megablast -query "$query" -db "$db" -evalue 1e-10 -max_target_seqs 200 -dust yes -outfmt "$OUTFMT" > "$out" 2>> "$log"
  fi
}

run_set(){
  local label=$1 kind=$2 fasta=$3 db=$4 final=$5
  local qdir="$S/tmp/${label}_queries" rdir="$S/tmp/${label}_results"
  split_fasta "$fasta" "$qdir" >/dev/null
  rm -rf "$rdir"; mkdir -p "$rdir"
  local q stem out log
  for q in "$qdir"/*.fa; do
    stem=$(basename "$q" .fa)
    out="$rdir/${stem}.tsv"; log="$S/logs/${label}_${stem}.log"
    if run_one_remote "$kind" "$q" "$db" "$out" "$log"; then
      echo "[PASS_REMOTE] $stem"
    elif run_one_local_if_available "$kind" "$q" "$db" "$out" "$log"; then
      echo "[PASS_LOCAL_FALLBACK] $stem"
    else
      echo "[FATAL] all transport routes failed for $label $stem; technical-censored, not a biological no-hit" >&2
      exit 2
    fi
  done
  : > "$final"
  for out in "$rdir"/*.tsv; do cat "$out" >> "$final"; done
}

echo '[1/3] RT5 vs NCBI nr, one query/request'
run_set RT5 protein "$S/work/RT5.faa" nr "$S/outputs/F3B4E_RT5_NCBI_NR.tsv"
echo '[2/3] partner5 vs NCBI nr, one query/request'
run_set PARTNER5 protein "$S/work/PARTNER5.faa" nr "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv"
echo '[3/3] architecture5 vs NCBI nt, one query/request'
run_set ARCH5 nucleotide "$S/work/ARCH5.fna" nt "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv"

META="$S/outputs/F3B4E_GLOBAL_REMOTE_RUNTIME_R1.json"
python - "$META" "$LOCK" <<'PY'
import json,sys,datetime
m,l=sys.argv[1:]
json.dump({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'transport_lock':json.load(open(l))},open(m,'w'),indent=2);open(m,'a').write('\n')
PY
PYTHONPATH="$S/scripts" python "$S/scripts/summarize_global.py" \
  "$S/work/TARGETS5.tsv" \
  "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" \
  "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" \
  "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" \
  "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" \
  "$S/work/BASE_LOCK.json" "$META" | tee "$S/logs/summarize_global_r1.log"

echo PASS_STAGE05F3B4E_G1_R1_REMOTE_TRANSPORT_HOTFIX_COMPLETE
