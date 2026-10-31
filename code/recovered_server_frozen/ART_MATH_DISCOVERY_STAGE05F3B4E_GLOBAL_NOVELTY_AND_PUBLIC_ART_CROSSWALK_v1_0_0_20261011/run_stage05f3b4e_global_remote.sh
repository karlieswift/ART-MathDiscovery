#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$S"/{outputs,logs,work,tmp}
for x in blastp blastn python; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x. Install: conda install -n art_novelty -y -c bioconda blast" >&2; exit 2; }; done
for f in "$S/inputs/F3B4E_RT5.faa" "$S/inputs/F3B4E_PARTNER5.faa" "$S/inputs/F3B4E_ARCHITECTURE5.fna" "$S/outputs/F3B4E_TARGETS_5.tsv" "$S/frozen/F3B4E_IMPLEMENTATION_LOCK.json"; do [[ -s "$f" ]] || { echo "[FATAL] run prepare first; missing $f" >&2; exit 2; }; done
OUTFMT='6 qseqid saccver pident length qstart qend sstart send evalue bitscore qlen slen qcovs stitle sscinames'
META="$S/outputs/F3B4E_GLOBAL_REMOTE_RUNTIME.json"
python - "$META" <<'PY'
import json,sys,datetime,subprocess
p=sys.argv[1]
def cmd(x):
 try:return subprocess.check_output(x,shell=True,text=True,stderr=subprocess.STDOUT).strip()
 except:return ''
json.dump({'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'blastp_version':cmd('blastp -version | head -1'),'blastn_version':cmd('blastn -version | head -1'),'databases':['NCBI nr remote','NCBI nt remote']},open(p,'w'),indent=2);open(p,'a').write('\n')
PY
run_retry(){
  local name=$1; shift
  local out="$S/outputs/$name" log="$S/logs/${name%.tsv}.log"
  for attempt in 1 2; do
    echo "[REMOTE] $name attempt=$attempt" | tee -a "$log"
    if "$@" > "$out" 2>> "$log"; then return 0; fi
    rm -f "$out"; sleep 20
  done
  echo "[FATAL] remote BLAST failed: $name" >&2; exit 2
}
echo '[1/3] RT5 vs NCBI nr remote BLASTP'
run_retry F3B4E_RT5_NCBI_NR.tsv blastp -query "$S/inputs/F3B4E_RT5.faa" -db nr -remote -evalue 1e-3 -max_target_seqs 200 -seg yes -outfmt "$OUTFMT"
echo '[2/3] partner5 vs NCBI nr remote BLASTP'
run_retry F3B4E_PARTNER5_NCBI_NR.tsv blastp -query "$S/inputs/F3B4E_PARTNER5.faa" -db nr -remote -evalue 1e-3 -max_target_seqs 200 -seg yes -outfmt "$OUTFMT"
echo '[3/3] architecture5 vs NCBI nt remote MEGABLAST'
run_retry F3B4E_ARCHITECTURE5_NCBI_NT.tsv blastn -task megablast -query "$S/inputs/F3B4E_ARCHITECTURE5.fna" -db nt -remote -evalue 1e-10 -max_target_seqs 200 -dust yes -outfmt "$OUTFMT"
python "$S/scripts/summarize_global.py" "$S/outputs/F3B4E_TARGETS_5.tsv" "$S/outputs/F3B4E_RT5_NCBI_NR.tsv" "$S/outputs/F3B4E_PARTNER5_NCBI_NR.tsv" "$S/outputs/F3B4E_ARCHITECTURE5_NCBI_NT.tsv" "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" "$S/frozen/F3B4E_IMPLEMENTATION_LOCK.json" "$META" | tee "$S/logs/summarize_global.log"
echo PASS_STAGE05F3B4E_GLOBAL_REMOTE_RUN_COMPLETE
