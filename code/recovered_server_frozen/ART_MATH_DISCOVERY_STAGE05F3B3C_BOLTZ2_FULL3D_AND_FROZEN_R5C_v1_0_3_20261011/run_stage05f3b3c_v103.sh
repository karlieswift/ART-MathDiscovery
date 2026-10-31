#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
S0="$P/ART_MATH_DISCOVERY_STAGE05F3B3C_BOLTZ2_FULL3D_AND_FROZEN_R5C_v1_0_2_20261011"
B="$P/ART_MATH_DISCOVERY_STAGE05F3B3B_KNOWN_SYSTEM_AUDIT_AND_STRUCTURE_HANDOFF_v1_0_0_20261011"
R5C="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007"
STATUS0="$S0/outputs/F3B3C_BOLTZ_STATUS.tsv"
AUD="$B/outputs/F3B3B_KNOWN_SYSTEM_AUDIT_32.tsv"
RTDB="$R5C/references/clean_rt_afsp"
LOCK="$S/frozen/F3B3C_V103_IMPLEMENTATION_LOCK.json"
mkdir -p "$S"/{frozen,raw/foldseek,outputs,logs,tmp/foldseek}
for x in "$STATUS0" "$AUD" "$S0/frozen/F3B3C_IMPLEMENTATION_LOCK.json" "$S/scripts/prepare_foldseek_jobs.py" "$S/scripts/summarize_r5c_v103.py"; do
  [[ -s "$x" ]] || { echo "[FATAL] missing $x" >&2; exit 2; }
done
conda run -n art_novelty bash -c 'command -v foldseek >/dev/null' >/dev/null 2>&1 || { echo '[FATAL] foldseek unavailable in art_novelty' >&2; exit 2; }
ls "${RTDB}"* >/dev/null 2>&1 || { echo '[FATAL] missing frozen R5C DB' >&2; exit 2; }

if [[ ! -s "$LOCK" ]]; then
  echo '[1/5] freeze v1.0.3 CRLF/Foldseek-handoff hotfix before Foldseek results'
  python - "$STATUS0" "$S0/frozen/F3B3C_IMPLEMENTATION_LOCK.json" "$S/scripts/prepare_foldseek_jobs.py" "$S/scripts/summarize_r5c_v103.py" "$S/run_stage05f3b3c_v103.sh" "$LOCK" <<'PY'
import csv,hashlib,json,sys
status,oldlock,p1,p2,runner,out=sys.argv[1:7]
def h(p):
    z=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):z.update(b)
    return z.hexdigest()
with open(status,newline='') as f:r=list(csv.DictReader(f,delimiter='\t'))
assert len(r)==32
bp=sum(x['boltz_status'].strip().startswith('PASS') for x in r)
bf=32-bp
L=json.load(open(oldlock));assert L['gate_pass'] is True
x={
 'stage':'Stage05F3B3C-v1.0.3-implementation-freeze',
 'gate':'PASS_F3B3C_V103_CRLF_FOLDSEEK_HANDOFF_HOTFIX_FROZEN_PRE_FOLDSEEK',
 'gate_pass':True,
 'hotfix_type':'IMPLEMENTATION_TEXT_SERIALIZATION_AND_FOLDSEEK_HANDOFF_ONLY',
 'v102_failure':'Python csv.DictWriter emitted CRLF; Bash read preserved carriage return in final cif_path field, causing first Foldseek query path to be invalid and the set -e pipeline to exit before any .m8 result',
 'v102_foldseek_result_n_before_hotfix':0,
 'boltz_status_sha256':h(status),
 'boltz_success_n_observed':bp,
 'boltz_technical_failure_n_observed':bf,
 'boltz_results_reused':True,
 'boltz_recomputed':False,
 'scientific_protocol_changed':False,
 'shortlist_membership_changed':False,
 'scientific_threshold_changed':False,
 'r5c_rules_changed':False,
 'prepare_jobs_sha256':h(p1),'summarizer_sha256':h(p2),'runner_sha256':h(runner),
 'claim_boundary':'Implementation handoff repair only; existing v1.0.2 Boltz results are reused byte-for-byte and frozen R5C rules are unchanged.'
}
with open(out,'w') as f:json.dump(x,f,indent=2);f.write('\n')
print(json.dumps(x,indent=2))
PY
  sha256sum "$LOCK" "$S/run_stage05f3b3c_v103.sh" "$S/scripts/prepare_foldseek_jobs.py" "$S/scripts/summarize_r5c_v103.py" > "$S/frozen/F3B3C_V103_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$S/run_stage05f3b3c_v103.sh" "$S/scripts/prepare_foldseek_jobs.py" "$S/scripts/summarize_r5c_v103.py"
else
  echo '[1/5] [SKIP] v1.0.3 implementation lock exists'
fi

echo '[2/5] prepare LF-normalized Foldseek jobs from frozen v1.0.2 Boltz status'
python "$S/scripts/prepare_foldseek_jobs.py" "$STATUS0" "$S/frozen/F3B3C_FOLDSEEK_JOBS.tsv" | tee "$S/logs/prepare_foldseek_jobs.log"

echo '[3/5] frozen R5C Foldseek searches via art_novelty'
printf 'candidate_id\tfoldseek_status\treturn_code\tm8_path\tlog_path\n' > "$S/outputs/F3B3C_FOLDSEEK_STATUS.tsv"
while IFS=$'\t' read -r rank fam cid acc len bst cif; do
  [[ "$rank" == shortlist_rank ]] && continue
  cif=${cif%$'\r'}
  out="$S/raw/foldseek/${cid}.m8"
  tmp="$S/tmp/foldseek/${cid}"
  log="$S/logs/foldseek_${cid}.log"
  rm -rf "$tmp"; rm -f "$out"
  set +e
  conda run -n art_novelty foldseek easy-search "$cif" "$RTDB" "$out" "$tmp" -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 --format-output 'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore' > "$log" 2>&1
  rc=$?
  set -e
  if [[ "$rc" -eq 0 && -e "$out" ]]; then
    st=PASS
    echo "[PASS] foldseek $cid"
  else
    st=TECHNICAL_FAILURE
    echo "[WARN] foldseek technical failure $cid rc=$rc" >&2
    tail -20 "$log" >&2 || true
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$cid" "$st" "$rc" "$out" "$log" >> "$S/outputs/F3B3C_FOLDSEEK_STATUS.tsv"
done < "$S/frozen/F3B3C_FOLDSEEK_JOBS.tsv"

echo '[4/5] frozen R5C classification with explicit Foldseek technical-failure handling'
python "$S/scripts/summarize_r5c_v103.py" "$STATUS0" "$S/outputs/F3B3C_FOLDSEEK_STATUS.tsv" "$S/raw/foldseek" "$AUD" "$S/outputs" "$LOCK" | tee "$S/logs/summarize_r5c.log"

echo '[5/5] hashes and compact result'
sha256sum "$LOCK" "$S/frozen/F3B3C_FOLDSEEK_JOBS.tsv" "$S/outputs/F3B3C_FOLDSEEK_STATUS.tsv" "$S/outputs/F3B3C_STRUCTURE_R5C_SUMMARY_32.tsv" "$S/outputs/F3B3C_ALL_CLEAN_RT_HITS.tsv" "$S/outputs/STAGE05F3B3C_GATE.json" > "$S/outputs/STAGE05F3B3C_V103_SHA256.tsv"
cat "$S/outputs/STAGE05F3B3C_GATE.json"
echo
echo '=== STRUCTURE/R5C SUMMARY 32 ==='
column -t -s $'\t' "$S/outputs/F3B3C_STRUCTURE_R5C_SUMMARY_32.tsv" | head -33 || true
echo
echo 'PASS_STAGE05F3B3C_V103_CRLF_FOLDSEEK_HANDOFF_HOTFIX_COMPLETE'
