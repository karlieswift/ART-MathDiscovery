#!/usr/bin/env bash
set -euo pipefail

THREADS=${1:-32}
NGPU=${2:-4}

S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
B="$P/ART_MATH_DISCOVERY_STAGE05F3B3B_KNOWN_SYSTEM_AUDIT_AND_STRUCTURE_HANDOFF_v1_0_0_20261011"
R5C="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007"

G="$B/outputs/STAGE05F3B3B_GATE.json"
QFA="$B/outputs/F3B3B_STRUCTURE_HANDOFF_32.faa"
AUD="$B/outputs/F3B3B_KNOWN_SYSTEM_AUDIT_32.tsv"
POL="$R5C/config/POLICY.json"
RTDB="$R5C/references/clean_rt_afsp"
LOCK="$S/frozen/F3B3C_IMPLEMENTATION_LOCK.json"

mkdir -p "$S"/{frozen,inputs,boltz_out,raw,outputs,logs,tmp}

# R1 implementation-only hotfix:
#   Boltz runs in art_structure.
#   Foldseek runs through the already-established art_novelty environment.
# Scientific universe, thresholds, shortlist and R5C rules are unchanged.
[[ "${CONDA_DEFAULT_ENV:-}" == "art_structure" ]] || {
  echo "[FATAL] activate art_structure before running this stage" >&2
  exit 2
}

for x in python sha256sum boltz conda; do
  command -v "$x" >/dev/null || { echo "[FATAL] missing $x in art_structure" >&2; exit 2; }
done

# Do not require Foldseek inside art_structure. Verify it in art_novelty.
conda run -n art_novelty bash -lc 'command -v foldseek >/dev/null' >/dev/null 2>&1 || {
  echo "[FATAL] foldseek not available in conda env art_novelty" >&2
  exit 2
}

run_foldseek() {
  conda run -n art_novelty foldseek "$@"
}

for x in \
  "$G" "$QFA" "$AUD" "$POL" \
  "$S/scripts/prepare_inputs.py" \
  "$S/scripts/run_boltz_workers.py" \
  "$S/scripts/summarize_r5c.py"
do
  [[ -s "$x" ]] || { echo "[FATAL] missing $x" >&2; exit 2; }
done

ls "${RTDB}"* >/dev/null 2>&1 || {
  echo "[FATAL] missing frozen R5C clean_rt_afsp DB" >&2
  exit 2
}

[[ "$NGPU" -ge 1 && "$NGPU" -le 4 ]] || {
  echo '[FATAL] NGPU must be 1..4' >&2
  exit 2
}

# Verify Boltz CLI in art_structure before freezing the implementation.
for flag in no_kernels use_msa_server seed recycling_steps sampling_steps diffusion_samples devices accelerator output_format out_dir; do
  boltz predict --help 2>&1 | grep -q -- "--${flag}" || {
    echo "[FATAL] installed boltz CLI lacks --${flag}" >&2
    exit 2
  }
done

if [[ ! -s "$LOCK" ]]; then
  echo '[1/7] freeze F3B3C v1.0.1 split-runtime implementation before structure results'

  find "$R5C/references" -maxdepth 1 -type f -name 'clean_rt_afsp*' -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "$S/frozen/R5C_CLEAN_RT_AFSP_SHA256.tsv"

  python - \
    "$G" "$QFA" "$AUD" "$POL" \
    "$S/frozen/R5C_CLEAN_RT_AFSP_SHA256.tsv" \
    "$S/scripts/prepare_inputs.py" \
    "$S/scripts/run_boltz_workers.py" \
    "$S/scripts/summarize_r5c.py" \
    "$S/run_stage05f3b3c.sh" \
    "$LOCK" <<'PY'
import hashlib, json, os, sys
(g,qfa,aud,pol,dbm,p1,p2,p3,runner,out)=sys.argv[1:11]

def h(p):
    x=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''):
            x.update(b)
    return x.hexdigest()

z=json.load(open(g))
assert z['gate_pass'] is True
assert z['structure_handoff_family_n']==32
assert z['shortlist_membership_changed'] is False
assert h(qfa)==z['shortlist_fasta_sha256']

lock={
  'stage':'Stage05F3B3C-v1.0.1-implementation-freeze',
  'gate':'PASS_F3B3C_V101_SPLIT_RUNTIME_HOTFIX_FROZEN_BEFORE_STRUCTURE_RESULTS',
  'gate_pass':True,
  'hotfix_type':'IMPLEMENTATION_ENVIRONMENT_SPLIT_ONLY',
  'v100_failure':'v1.0.0 required foldseek in art_structure during preflight; no structure computation or scientific result occurred before failure',
  'scientific_protocol_changed':False,
  'f3b3b_gate_sha256':h(g),
  'structure_handoff_fasta_sha256':h(qfa),
  'known_audit_sha256':h(aud),
  'r5c_policy_sha256':h(pol),
  'r5c_clean_rt_db_hash_manifest_sha256':h(dbm),
  'structure_handoff_family_n':32,
  'boltz_runtime_env':'art_structure',
  'foldseek_runtime_env':'art_novelty via conda run',
  'boltz_model':'boltz2',
  'msa_mode':'server-assisted',
  'seed':20261007,
  'recycling_steps':3,
  'sampling_steps':200,
  'diffusion_samples':1,
  'no_kernels':True,
  'output_format':'mmcif',
  'gpu_execution':'independent single-GPU workers; no multi-GPU DDP/NCCL',
  'max_gpu_workers':4,
  'no_single_sequence_fallback':True,
  'no_length_truncation':True,
  'technical_failure_policy':'record unresolved; do not change scientific settings or silently substitute msa: empty',
  'foldseek_search_evalue':1.0,
  'foldseek_max_seqs':100,
  'foldseek_exact_tmscore':1,
  'foldseek_format':'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore',
  'strong_rule':'alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50',
  'suggestive_rule':'alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30',
  'shortlist_membership_changed':False,
  'known_matches_removed':False,
  'novelty_used_for_structure_handoff':False,
  'scientific_threshold_changed':False,
  'prepare_inputs_sha256':h(p1),
  'boltz_worker_sha256':h(p2),
  'summarizer_sha256':h(p3),
  'runner_sha256':h(runner),
  'claim_boundary':'Frozen full-3D/clean-RT structural compatibility test only; no biochemical RT/ART or novelty claim is authorized from this stage alone.'
}

with open(out,'w') as f:
    json.dump(lock,f,indent=2)
    f.write('\n')
print(json.dumps(lock,indent=2))
PY

  sha256sum \
    "$LOCK" \
    "$S/run_stage05f3b3c.sh" \
    "$S/scripts/prepare_inputs.py" \
    "$S/scripts/run_boltz_workers.py" \
    "$S/scripts/summarize_r5c.py" \
    > "$S/frozen/F3B3C_V101_FREEZE_SHA256.tsv"

  chmod a-w \
    "$LOCK" \
    "$S/run_stage05f3b3c.sh" \
    "$S/scripts/prepare_inputs.py" \
    "$S/scripts/run_boltz_workers.py" \
    "$S/scripts/summarize_r5c.py"
else
  echo '[1/7] [SKIP] v1.0.1 implementation lock exists'
fi

python - \
  "$LOCK" \
  "$S/run_stage05f3b3c.sh" \
  "$S/scripts/prepare_inputs.py" \
  "$S/scripts/run_boltz_workers.py" \
  "$S/scripts/summarize_r5c.py" <<'PY'
import hashlib,json,sys
L=json.load(open(sys.argv[1]))

def h(p):
    x=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            x.update(b)
    return x.hexdigest()

assert h(sys.argv[2])==L['runner_sha256']
assert h(sys.argv[3])==L['prepare_inputs_sha256']
assert h(sys.argv[4])==L['boltz_worker_sha256']
assert h(sys.argv[5])==L['summarizer_sha256']
assert L['scientific_protocol_changed'] is False
print('[PASS] F3B3C v1.0.1 frozen split-runtime implementation verified')
PY

echo '[2/7] prepare exact frozen 32 Boltz inputs'
python \
  "$S/scripts/prepare_inputs.py" \
  "$QFA" "$AUD" "$S/inputs" \
  "$S/frozen/F3B3C_QUERY_MANIFEST.tsv" \
  | tee "$S/logs/prepare_inputs.log"

echo '[3/7] Boltz-2 MSA-assisted full-3D primary attempt in art_structure'
export BOLTZ_CACHE=${BOLTZ_CACHE:-"$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache"}
mkdir -p "$BOLTZ_CACHE"

python \
  "$S/scripts/run_boltz_workers.py" \
  "$S/frozen/F3B3C_QUERY_MANIFEST.tsv" \
  "$S/boltz_out" \
  "$S/logs" \
  "$S/outputs/F3B3C_BOLTZ_STATUS.tsv" \
  "$NGPU" \
  | tee "$S/logs/boltz_workers_summary.log"

echo '[4/7] exact frozen R5C Foldseek searches via art_novelty'
mkdir -p "$S/raw/foldseek" "$S/tmp/foldseek"

tail -n +2 "$S/outputs/F3B3C_BOLTZ_STATUS.tsv" \
| while IFS=$'\t' read -r rank fam cid acc len gpu bst rc elapsed cif; do
    [[ "$bst" == PASS* ]] || continue

    out="$S/raw/foldseek/${cid}.m8"
    tmp="$S/tmp/foldseek/${cid}"

    rm -rf "$tmp"
    rm -f "$out"

    run_foldseek easy-search \
      "$cif" \
      "$RTDB" \
      "$out" \
      "$tmp" \
      -e 1.0 \
      --max-seqs 100 \
      --threads "$THREADS" \
      --exact-tmscore 1 \
      --format-output 'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore' \
      > "$S/logs/foldseek_${cid}.log" 2>&1

    echo "[PASS] foldseek $cid"
  done

echo '[5/7] frozen R5C classification'
python \
  "$S/scripts/summarize_r5c.py" \
  "$S/outputs/F3B3C_BOLTZ_STATUS.tsv" \
  "$S/raw/foldseek" \
  "$AUD" \
  "$S/outputs" \
  "$LOCK" \
  | tee "$S/logs/summarize_r5c.log"

echo '[6/7] hashes'
sha256sum \
  "$LOCK" \
  "$S/frozen/F3B3C_QUERY_MANIFEST.tsv" \
  "$S/outputs/F3B3C_BOLTZ_STATUS.tsv" \
  "$S/outputs/F3B3C_STRUCTURE_R5C_SUMMARY_32.tsv" \
  "$S/outputs/F3B3C_ALL_CLEAN_RT_HITS.tsv" \
  "$S/outputs/STAGE05F3B3C_GATE.json" \
  > "$S/outputs/STAGE05F3B3C_SHA256.tsv"

echo '[7/7] compact result'
cat "$S/outputs/STAGE05F3B3C_GATE.json"
echo
echo '=== STRUCTURE/R5C SUMMARY 32 ==='
column -t -s $'\t' "$S/outputs/F3B3C_STRUCTURE_R5C_SUMMARY_32.tsv" | head -33 || true
echo
echo 'PASS_STAGE05F3B3C_V101_SPLIT_RUNTIME_PRIMARY_STRUCTURE_ATTEMPT_COMPLETE'
