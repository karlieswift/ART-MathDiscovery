#!/usr/bin/env bash
set -euo pipefail

S=$(cd "$(dirname "$0")/.." && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1

SHARD=${1:?usage: run_stage03_shard_r2h2.sh SHARD.txt [CPU] [JOBS]}
CPU=${2:-12}
JOBS=${3:-4}

EXPECTED=74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e

HMM="${RVT_HMM:-$P/ART_MATH_DISCOVERY_STAGE05D1_NOVELTY_FIRST_MINING_v1_0_0_20261006/references/myRT/RVT-All.hmm}"

for x in zstd hmmscan python; do
    command -v "$x" >/dev/null || {
        echo "[FATAL] missing $x" >&2
        exit 2
    }
done

[[ -s "$SHARD" ]] || {
    echo "[FATAL] missing shard $SHARD" >&2
    exit 2
}

[[ -s "$HMM" ]] || {
    echo "[FATAL] missing RVT-All.hmm" >&2
    exit 2
}

ACTUAL=$(sha256sum "$HMM" | awk '{print $1}')
[[ "$ACTUAL" == "$EXPECTED" ]] || {
    echo "[FATAL] RVT-All.hmm SHA mismatch: $ACTUAL" >&2
    exit 2
}

[[ -x "$S/scripts/scan_one_accession.sh" ]] || exit 2
[[ -s "$S/scripts/extract_fasta_ids.py" ]] || exit 2

mkdir -p \
  "$S/work/per_accession" \
  "$S/work/done" \
  "$S/outputs/formal_faa" \
  "$S/logs/per_accession"

export KEEP_ORF_FASTA=0

echo "[R2H2 SHARD START]"
echo "host=$(hostname)"
echo "shard=$SHARD"
echo "accession_n=$(wc -l < "$SHARD")"
echo "cpu_per_hmmscan=$CPU"
echo "parallel_jobs=$JOBS"
echo "rvt_hmm_sha256=$ACTUAL"

cat "$SHARD" |
  xargs -P "$JOBS" -I{} \
    bash "$S/scripts/scan_one_accession.sh" \
      "$S" "{}" "$HMM" "$CPU"

echo "PASS_STAGE05F3V2_STAGE03_SHARD_COMPLETE host=$(hostname) shard=$(basename "$SHARD")"
