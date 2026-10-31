#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
NGPU=${2:-4}
S=$(cd "$(dirname "$0")" && pwd)
P=${3:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
B2="$P/ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B2_HISTORICAL_F3B3A_EXACT_PRIORITY_PROJECTION_v1_0_1_20261016"
HIST="$P/ART_MATH_DISCOVERY_STAGE05F3B3C_BOLTZ2_FULL3D_AND_FROZEN_R5C_v1_0_2_20261011"
R5C="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007"
G="$B2/outputs/STAGE04B2_GATE.json"
QTSV="$B2/outputs/F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.tsv"
QFA="$B2/outputs/F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.faa"
POL="$R5C/config/POLICY.json"
RTDB="$R5C/references/clean_rt_afsp"
HIST_LOCK="$HIST/frozen/F3B3C_IMPLEMENTATION_LOCK.json"
HIST_GATE="$HIST/outputs/STAGE05F3B3C_GATE.json"
LOCK="$S/frozen/STAGE04B3_IMPLEMENTATION_AND_INPUT_LOCK.json"
HFA="$S/inputs/F3V2_STAGE04B3_STRUCTURE_HANDOFF_32.faa"
HAUD="$S/inputs/F3V2_STAGE04B3_HANDOFF_AUDIT_32.tsv"
HOUT="$S/outputs/historical_r5c_exact"
mkdir -p "$S"/{frozen,inputs,boltz_out,raw/foldseek,tmp/foldseek,outputs,logs} "$HOUT"

echo '[0/8] selftests and runtime preflight'
python "$S/scripts/selftest.py"
python -m py_compile "$S"/scripts/*.py
echo PASS_STAGE04B3_PYTHON_COMPILE_SELFTEST
for x in python sha256sum boltz conda; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
python - <<'PY'
import torch,sys
print('[RUNTIME] python_prefix =',sys.prefix)
print('[RUNTIME] torch =',torch.__version__)
print('[RUNTIME] cuda_available =',torch.cuda.is_available())
print('[RUNTIME] cuda_device_count =',torch.cuda.device_count())
if not torch.cuda.is_available(): raise SystemExit('[FATAL] CUDA unavailable')
PY
conda run -n art_novelty bash -c 'command -v foldseek >/dev/null' >/dev/null 2>&1 || { echo '[FATAL] foldseek missing in art_novelty' >&2; exit 2; }
[[ "$NGPU" -ge 1 && "$NGPU" -le 4 ]] || { echo '[FATAL] NGPU must be 1..4' >&2; exit 2; }
for x in "$G" "$QTSV" "$QFA" "$POL" "$HIST_LOCK"; do [[ -s "$x" ]] || { echo "[FATAL] missing $x" >&2; exit 2; }; done
ls "${RTDB}"* >/dev/null 2>&1 || { echo '[FATAL] clean_rt_afsp DB missing' >&2; exit 2; }
BOLTZ_HELP=$(boltz predict --help 2>&1)
for flag in no_kernels use_msa_server seed recycling_steps sampling_steps diffusion_samples devices accelerator output_format out_dir; do grep -q -- "--${flag}" <<<"$BOLTZ_HELP" || { echo "[FATAL] boltz CLI lacks --$flag" >&2; exit 2; }; done

echo '[1/8] exact Stage04B2 delta32 handoff'
python "$S/scripts/prepare_delta32_handoff.py" "$G" "$QTSV" "$QFA" "$HFA" "$HAUD" | tee "$S/logs/prepare_handoff.log"

echo '[2/8] freeze exact historical structure/R5C settings before delta structure results'
find "$R5C/references" -maxdepth 1 -type f -name 'clean_rt_afsp*' -print0 | sort -z | xargs -0 sha256sum > "$S/frozen/R5C_CLEAN_RT_AFSP_SHA256.tsv"
python - "$G" "$QTSV" "$QFA" "$HFA" "$HAUD" "$POL" "$HIST_LOCK" "$S/frozen/R5C_CLEAN_RT_AFSP_SHA256.tsv" "$S/scripts/prepare_inputs_historical_exact.py" "$S/scripts/run_boltz_workers_historical_exact.py" "$S/scripts/summarize_r5c_historical_exact.py" "$S/scripts/prepare_delta32_handoff.py" "$S/run_stage04b3_delta32.sh" "$LOCK" <<'PY'
import hashlib,json,sys
(g,qtsv,qfa,hfa,haud,pol,hlock,dbm,p1,p2,p3,adapter,runner,out)=sys.argv[1:15]
def h(p):
 x=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(16<<20),b''): x.update(b)
 return x.hexdigest()
b2=json.load(open(g)); H=json.load(open(hlock))
assert b2['gate_pass'] is True and b2['structure_priority_shortlist_family_n']==32
assert H['gate_pass'] is True and H['structure_handoff_family_n']==32
assert H['boltz_model']=='boltz2' and H['seed']==20261007 and H['recycling_steps']==3 and H['sampling_steps']==200 and H['diffusion_samples']==1 and H['no_kernels'] is True
assert H['foldseek_search_evalue']==1.0 and H['foldseek_max_seqs']==100 and H['foldseek_exact_tmscore']==1
assert H['strong_rule']=='alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50'
assert H['suggestive_rule']=='alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30'
assert H['no_single_sequence_fallback'] is True and H['no_length_truncation'] is True
assert h(p1)=='71845ec178195cf585c6ebb27ab0d1c85e261e31f3358db705240ca35d86d901'
assert h(p2)=='8d53fae473fdba5724f74cf1b5ca422eb2bec08bd9654009c66adc227949f4af'
assert h(p3)=='27418c798226f7a61464b9f7f3e12e264bd3b33320247187baf5c1327facf93b'
L={
 'stage':'Stage05F3V2-Stage04B3-delta32-implementation-freeze','gate':'PASS_F3V2_STAGE04B3_EXACT_HISTORICAL_F3B3C_PROTOCOL_FROZEN_BEFORE_DELTA_STRUCTURE_RESULTS','gate_pass':True,
 'historical_source_stage':'Stage05F3B3C v1.0.2 (2026-10-11)','historical_structure_archive_sha256':'0ebc19cd4628fffcce0a9fc67de7fdd9b80571b44bf291f0ddd5327868c40dc5',
 'historical_runtime_lock_sha256':h(hlock),'stage04b2_gate_sha256':h(g),'stage04b2_shortlist_tsv_sha256':h(qtsv),'stage04b2_shortlist_fasta_sha256':h(qfa),
 'delta32_handoff_fasta_sha256':h(hfa),'delta32_handoff_audit_sha256':h(haud),'r5c_policy_sha256':h(pol),'r5c_clean_rt_db_hash_manifest_sha256':h(dbm),
 'historical_prepare_inputs_sha256':h(p1),'historical_boltz_worker_sha256':h(p2),'historical_summarizer_sha256':h(p3),'delta_handoff_adapter_sha256':h(adapter),'runner_sha256':h(runner),
 'structure_handoff_family_n':32,'boltz_model':'boltz2','msa_mode':'server-assisted','seed':20261007,'recycling_steps':3,'sampling_steps':200,'diffusion_samples':1,'no_kernels':True,'output_format':'mmcif',
 'gpu_execution':'independent single-GPU workers; no multi-GPU DDP/NCCL','max_gpu_workers':4,'no_single_sequence_fallback':True,'no_length_truncation':True,
 'technical_failure_policy':'record unresolved; do not change scientific settings or silently substitute msa: empty','foldseek_search_evalue':1.0,'foldseek_max_seqs':100,'foldseek_exact_tmscore':1,
 'foldseek_format':'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore','strong_rule':'alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50','suggestive_rule':'alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30',
 'shortlist_membership_changed':False,'known_sequence_audit_used_for_handoff':False,'novelty_used_for_handoff':False,'scientific_threshold_changed':False,'science_change':False,
 'claim_boundary':'Out-of-universe structure generalization only. Exact Stage04B2 delta32 membership is tested with the exact historical F3B3C Boltz-2/R5C protocol; no membership retuning or known-system/novelty filtering is authorized.'}
with open(out,'w') as f: json.dump(L,f,indent=2); f.write('\n')
print(json.dumps(L,indent=2))
PY

echo '[3/8] prepare exact historical Boltz input schema'
python "$S/scripts/prepare_inputs_historical_exact.py" "$HFA" "$HAUD" "$S/inputs/boltz_yaml" "$S/frozen/F3V2_STAGE04B3_QUERY_MANIFEST.tsv" | tee "$S/logs/prepare_inputs.log"

echo '[4/8] Boltz-2 MSA-assisted full-3D on 32 frozen delta representatives'
export BOLTZ_CACHE=${BOLTZ_CACHE:-"$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache"}
mkdir -p "$BOLTZ_CACHE"
python "$S/scripts/run_boltz_workers_historical_exact.py" "$S/frozen/F3V2_STAGE04B3_QUERY_MANIFEST.tsv" "$S/boltz_out" "$S/logs" "$S/outputs/F3V2_STAGE04B3_BOLTZ_STATUS.tsv" "$NGPU" | tee "$S/logs/boltz_workers_summary.log"

echo '[5/8] exact historical R5C Foldseek searches via art_novelty'
run_foldseek(){ conda run -n art_novelty foldseek "$@"; }
tail -n +2 "$S/outputs/F3V2_STAGE04B3_BOLTZ_STATUS.tsv" | while IFS=$'\t' read -r rank fam cid acc len gpu bst rc elapsed cif; do
 [[ "$bst" == PASS* ]] || continue
 out="$S/raw/foldseek/${cid}.m8"; tmp="$S/tmp/foldseek/${cid}"
 rm -rf "$tmp"; rm -f "$out"
 run_foldseek easy-search "$cif" "$RTDB" "$out" "$tmp" -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 --format-output 'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore' > "$S/logs/foldseek_${cid}.log" 2>&1
 echo "[PASS] foldseek $cid"
done

echo '[6/8] exact historical R5C classification'
python "$S/scripts/summarize_r5c_historical_exact.py" "$S/outputs/F3V2_STAGE04B3_BOLTZ_STATUS.tsv" "$S/raw/foldseek" "$HAUD" "$HOUT" "$LOCK" | tee "$S/logs/summarize_r5c.log"

echo '[7/8] Stage04B3 delta-generalization gate'
python "$S/scripts/finalize_delta32.py" "$HOUT/STAGE05F3B3C_GATE.json" "$HOUT/F3B3C_STRUCTURE_R5C_SUMMARY_32.tsv" "$HOUT/F3B3C_ALL_CLEAN_RT_HITS.tsv" "$G" "$LOCK" "$S/outputs" "$HIST_GATE" | tee "$S/logs/finalize_delta32.log"

echo '[8/8] hashes and compact result'
sha256sum "$LOCK" "$S/frozen/F3V2_STAGE04B3_QUERY_MANIFEST.tsv" "$S/outputs/F3V2_STAGE04B3_BOLTZ_STATUS.tsv" "$S/outputs/F3V2_STAGE04B3_STRUCTURE_R5C_SUMMARY_32.tsv" "$S/outputs/F3V2_STAGE04B3_ALL_CLEAN_RT_HITS.tsv" "$S/outputs/STAGE04B3_GATE.json" > "$S/outputs/STAGE04B3_SHA256.tsv"
cat "$S/outputs/STAGE04B3_GATE.json"
echo
echo '=== DELTA32 STRUCTURE/R5C SUMMARY ==='
column -t -s $'\t' "$S/outputs/F3V2_STAGE04B3_STRUCTURE_R5C_SUMMARY_32.tsv" | head -33 || true
echo
if python - "$S/outputs/STAGE04B3_GATE.json" <<'PY'
import json,sys
z=json.load(open(sys.argv[1])); raise SystemExit(0 if z['gate_pass'] else 1)
PY
then echo PASS_STAGE05F3V2_STAGE04B3_DELTA32_FROZEN_STRUCTURE_GENERALIZATION_PRIMARY_COMPLETE; else exit 2; fi
