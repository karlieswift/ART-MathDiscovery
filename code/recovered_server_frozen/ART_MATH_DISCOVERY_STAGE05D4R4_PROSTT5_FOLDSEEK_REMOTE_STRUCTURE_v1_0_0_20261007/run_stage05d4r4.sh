#!/usr/bin/env bash
set -euo pipefail
R3=${1:?usage: run_stage05d4r4.sh STAGE05D4R3_ROOT [THREADS] [GPU]}
THREADS=${2:-32}
GPU=${3:-0}
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
mkdir -p outputs raw tmp/search logs
python - "$R3" <<'PY'
import json,sys,os
r=sys.argv[1]
g=json.load(open(os.path.join(r,'outputs','GATE_STAGE05D4R3.json')))
assert g['gate_pass'] and g['gate']=='PASS_COMPACT_STRUCTURE_QUERY_TIERS_FROZEN',g
assert g['n_tier_A']==350,g
print('[PASS] R3 gate and Tier-A count')
PY
command -v foldseek >/dev/null || { echo '[ERROR] foldseek not found'; exit 2; }
[ -e references/prostt5 ] || { echo '[ERROR] references/prostt5 missing; run ./prepare_foldseek_refs.sh'; exit 2; }
[ -e references/af_swissprot ] || { echo '[ERROR] references/af_swissprot missing; run ./prepare_foldseek_refs.sh'; exit 2; }
cp -f "$R3/outputs/STAGE05D4R3_TIER_A_COMPACT_STRUCTURE_READY.faa" work/TIER_A.faa
cp -f "$R3/outputs/STAGE05D4R3_NR90_STRUCTURE_PRIORITY.tsv" work/R3_PRIORITY.tsv
ARGS=(easy-search work/TIER_A.faa references/af_swissprot raw/TIER_A_AFSP.m8 tmp/search --prostt5-model references/prostt5/prostt5-f16.gguf -e 0.01 --max-seqs 50 --threads "$THREADS" --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident)
if [ "$GPU" = "1" ]; then ARGS+=(--gpu 1); fi
printf '[RUN] foldseek %q ' "${ARGS[@]}"; echo
foldseek "${ARGS[@]}" 2>&1 | tee logs/foldseek.log
SP=NA
[ -s references/uniprot_sprot.fasta.gz ] && SP=references/uniprot_sprot.fasta.gz
python scripts/summarize_foldseek.py raw/TIER_A_AFSP.m8 work/R3_PRIORITY.tsv work/TIER_A.faa outputs "$SP"
python scripts/write_gate.py outputs
sha256sum raw/TIER_A_AFSP.m8 outputs/STAGE05D4R4_TIER_A_STRUCTURAL_TRIAGE.tsv outputs/STAGE05D4R4_SUMMARY.json outputs/GATE_STAGE05D4R4.json > outputs/STAGE05D4R4_SHA256.tsv
