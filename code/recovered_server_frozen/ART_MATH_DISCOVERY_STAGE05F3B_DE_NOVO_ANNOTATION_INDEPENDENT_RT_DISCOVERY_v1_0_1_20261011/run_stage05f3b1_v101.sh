#!/usr/bin/env bash
set -euo pipefail

S=$(cd "$(dirname "$0")" && pwd)

LOCK="$S/frozen/F3B0_V101_IMPLEMENTATION_LOCK.json"
ORF="$S/work/F3B1_ALL_ORFS_NUMERIC.faa"
HMM="$S/references/RVT-All.hmm"

ORF_N=65315646

PARALLEL=$(
python - "$LOCK" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))["parallel_processes"])
PY
)

SHARDS=$(
python - "$LOCK" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))["shard_n"])
PY
)

python - "$LOCK" "$S/run_stage05f3b1_v101.sh" \
  "$S/scripts/shard_fasta.py" "$S/scripts/collect_f3b1.py" <<'PY'
import hashlib,json,sys

lock,runner,sharder,collector=sys.argv[1:5]
x=json.load(open(lock))

def h(p):
    q=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            q.update(b)
    return q.hexdigest()

assert h(runner)==x["runner_sha256"]
assert h(sharder)==x["sharder_sha256"]
assert h(collector)==x["collector_sha256"]

print("[PASS] F3B0 v1.0.1 implementation lock verified")
PY

echo "[1/4] shard 65,315,646 frozen ORFs"

if [[ ! -s "$S/frozen/SHARD_SPLIT_COMPLETE.json" ]]; then
    python "$S/scripts/shard_fasta.py" \
      "$ORF" \
      "$S/work/shards" \
      "$SHARDS" \
      "$ORF_N" \
      "$S/frozen/SHARD_COUNTS.tsv" \
      "$S/frozen/SHARD_SPLIT_COMPLETE.json"
else
    echo "[SKIP] frozen shards already complete"
fi

SUM=$(
  awk 'NR>1{s+=$2}END{print s+0}' \
    "$S/frozen/SHARD_COUNTS.tsv"
)

[[ "$SUM" -eq "$ORF_N" ]] || {
    echo "[FATAL] shard count mismatch: $SUM"
    exit 2
}

echo "[2/4] parallel frozen RVT-All hmmscan"
echo "parallel_processes=$PARALLEL shards=$SHARDS"

export S HMM

cut -f1 "$S/frozen/SHARD_COUNTS.tsv" |
tail -n +2 |
xargs -P "$PARALLEL" -I{} \
bash -c '
set -euo pipefail

SH="$1"
S="$2"
HMM="$3"

DONE="$S/raw/shards/${SH}.done"

if [[ -s "$DONE" ]]; then
    echo "[SKIP] $SH"
    exit 0
fi

TBL="$S/raw/shards/${SH}.tblout"
DOM="$S/raw/shards/${SH}.domtblout"

TT="${TBL}.part"
DD="${DOM}.part"

rm -f "$TT" "$DD"

OMP_NUM_THREADS=1 \
OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 \
hmmscan \
  --cpu 1 \
  --noali \
  -o /dev/null \
  -E 1e-2 \
  --domE 1000 \
  --tblout "$TT" \
  --domtblout "$DD" \
  "$HMM" \
  "$S/work/shards/${SH}.faa" \
  > "$S/logs/shards/${SH}.stdout" \
  2> "$S/logs/shards/${SH}.stderr"

mv "$TT" "$TBL"
mv "$DD" "$DOM"

printf "PASS\n" > "$DONE"

echo "[PASS] $SH"
' _ {} "$S" "$HMM"

DONE_N=$(
  find "$S/raw/shards" \
    -maxdepth 1 \
    -name 'shard_*.done' \
    -type f | wc -l
)

[[ "$DONE_N" -eq "$SHARDS" ]] || {
    echo "[FATAL] completed shards=$DONE_N expected=$SHARDS"
    exit 2
}

echo "[3/4] merge shard HMM outputs"

MERGED_T="$S/raw/STAGE05F3B1_RVT_ALL.tblout"
MERGED_D="$S/raw/STAGE05F3B1_RVT_ALL.domtblout"

: > "$MERGED_T"
: > "$MERGED_D"

while IFS=$'\t' read -r SH N B
do
    [[ "$SH" == "shard" ]] && continue

    cat "$S/raw/shards/${SH}.tblout" >> "$MERGED_T"
    cat "$S/raw/shards/${SH}.domtblout" >> "$MERGED_D"

done < "$S/frozen/SHARD_COUNTS.tsv"

echo "[4/4] frozen classification"

python \
  "$S/scripts/collect_f3b1.py" \
  "$MERGED_T" \
  "$MERGED_D" \
  "$LOCK" \
  "$S/outputs" \
  "$ORF_N" \
  | tee "$S/logs/collect_f3b1.log"

sha256sum \
  "$LOCK" \
  "$MERGED_T" \
  "$MERGED_D" \
  "$S/outputs/STAGE05F3B1_GATE.json" \
  "$S/outputs/STAGE05F3B1_CANDIDATE_SUMMARY.tsv" \
  "$S/outputs/STAGE05F3B1_RVT_PAIR_HITS.tsv.gz" \
  > "$S/outputs/STAGE05F3B1_V101_SHA256.tsv"

echo
echo "PASS_STAGE05F3B1_V101_PARALLEL_EQUIVALENT_SCREEN_COMPLETE"
