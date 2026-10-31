#!/usr/bin/env bash
set -euo pipefail

THREADS=${1:-32}

S=$(cd "$(dirname "$0")" && pwd)

P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1

SRC="$P/external_db/F3_LARGE_UNIVERSE_20261009/F3A_EXTERNAL_UNIVERSE.fna"

LOCK="$S/frozen/F3B0_IMPLEMENTATION_LOCK.json"

HMM="$S/references/RVT-All.hmm"

EXPECTED_SOURCE=acea0963b7f60eb6d273808c2ef1d73d833a654870faa301752e491a3759ef54

mkdir -p "$S"/{work,raw,outputs,logs,tmp}

for x in mmseqs hmmscan sha256sum python; do
    command -v "$x" >/dev/null || {
        echo "[FATAL] missing $x"
        exit 2
    }
done

ACTUAL=$(sha256sum "$SRC" | awk '{print $1}')

[[ "$ACTUAL" == "$EXPECTED_SOURCE" ]] || {
    echo "[FATAL] F3B source lock mismatch"
    exit 2
}

python - "$LOCK" "$S/run_stage05f3b1.sh" "$S/scripts/collect_f3b1.py" <<'PY'
import json,hashlib,sys

lock,runner,collector=sys.argv[1:4]
x=json.load(open(lock))

def h(p):
    q=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            q.update(b)
    return q.hexdigest()

assert h(runner)==x["runner_sha256"]
assert h(collector)==x["collector_sha256"]

print("[PASS] F3B0 implementation lock verified")
PY

echo "[1/5] nucleotide DB"

rm -rf \
  "$S/work/genomeDB"* \
  "$S/work/F3B1_orfs_aa"* \
  "$S/work/F3B1_ALL_ORFS_NUMERIC.faa"

mmseqs createdb \
  "$SRC" \
  "$S/work/genomeDB" \
  --shuffle 0 \
  > "$S/logs/mmseqs_createdb.log" 2>&1

echo "[2/5] six-frame maximal stop-to-stop translation"

mmseqs extractorfs \
  "$S/work/genomeDB" \
  "$S/work/F3B1_orfs_aa" \
  --orf-start-mode 1 \
  --forward-frames 1,2,3 \
  --reverse-frames 1,2,3 \
  --contig-start-mode 2 \
  --contig-end-mode 2 \
  --min-length 45 \
  --max-length 1000000 \
  --max-gaps 2147483647 \
  --translation-table 11 \
  --translate 1 \
  --create-lookup 1 \
  --threads "$THREADS" \
  > "$S/logs/mmseqs_extractorfs.log" 2>&1

ORF_N=$(wc -l < "$S/work/F3B1_orfs_aa.index")

echo "F3B1_ORF_N=$ORF_N"

echo "[3/5] ORF FASTA"

mmseqs convert2fasta \
  "$S/work/F3B1_orfs_aa" \
  "$S/work/F3B1_ALL_ORFS_NUMERIC.faa" \
  --use-header-file 0 \
  > "$S/logs/mmseqs_convert2fasta.log" 2>&1

echo "[4/5] frozen RVT-All hmmscan"

hmmscan \
  --cpu "$THREADS" \
  --noali \
  -E 1e-2 \
  --domE 1000 \
  --tblout "$S/raw/STAGE05F3B1_RVT_ALL.tblout" \
  --domtblout "$S/raw/STAGE05F3B1_RVT_ALL.domtblout" \
  "$HMM" \
  "$S/work/F3B1_ALL_ORFS_NUMERIC.faa" \
  > "$S/logs/hmmscan_RVT_ALL.stdout" \
  2> "$S/logs/hmmscan_RVT_ALL.stderr"

echo "[5/5] frozen classification"

python \
  "$S/scripts/collect_f3b1.py" \
  "$S/raw/STAGE05F3B1_RVT_ALL.tblout" \
  "$S/raw/STAGE05F3B1_RVT_ALL.domtblout" \
  "$LOCK" \
  "$S/outputs" \
  "$ORF_N" \
  | tee "$S/logs/collect_f3b1.log"

if [[ -s "$S/outputs/STAGE05F3B1_FORMAL_IDS.txt" ]]; then
    mmseqs createsubdb \
      "$S/outputs/STAGE05F3B1_FORMAL_IDS.txt" \
      "$S/work/F3B1_orfs_aa" \
      "$S/outputs/F3B1_FORMAL_DB"

    mmseqs convert2fasta \
      "$S/outputs/F3B1_FORMAL_DB" \
      "$S/outputs/STAGE05F3B1_FORMAL_CANDIDATES.faa" \
      --use-header-file 1
fi

sha256sum \
  "$S/frozen/F3B0_IMPLEMENTATION_LOCK.json" \
  "$S/outputs/STAGE05F3B1_GATE.json" \
  "$S/outputs/STAGE05F3B1_CANDIDATE_SUMMARY.tsv" \
  "$S/outputs/STAGE05F3B1_RVT_PAIR_HITS.tsv.gz" \
  > "$S/outputs/STAGE05F3B1_SHA256.tsv"

echo
echo "[PASS] Stage05F3B1 complete"
