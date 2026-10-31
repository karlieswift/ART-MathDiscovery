#!/usr/bin/env bash
set -euo pipefail

P=${1:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
THREADS=${2:-64}

OLD="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_1_20261011"
DELTA="$P/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_DISTRIBUTED_ACQUISITION_R2_v1_0_2_20261013"
HERE=$(cd "$(dirname "$0")" && pwd)
OUT="$HERE/outputs"
WORK="$HERE/work"
TMP="$WORK/mmseqs_tmp"

mkdir -p "$OUT" "$WORK" "$TMP"

command -v python >/dev/null
command -v mmseqs >/dev/null

MMSEQS=$(command -v mmseqs)
MMSEQS_VERSION=$(mmseqs version | head -n1)

printf '%s\n' "[STAGE04A START]" \
  "project_root=$P" \
  "old_root=$OLD" \
  "delta_root=$DELTA" \
  "threads=$THREADS" \
  "mmseqs=$MMSEQS" \
  "mmseqs_version=$MMSEQS_VERSION"

python "$HERE/scripts/prepare_combined_exact_universe.py" \
  --project-root "$P" \
  --old-root "$OLD" \
  --delta-root "$DELTA" \
  --out "$OUT" \
  --old-expected 10329 \
  --delta-expected 85503

rm -rf "$TMP"/*
rm -f "$OUT"/combined60g_cluster.tsv "$OUT"/combined60g_rep_seq.fasta "$OUT"/combined60g_all_seqs.fasta

mmseqs easy-cluster \
  "$OUT/COMBINED_GLOBAL_EXACT_UNIQUE.faa" \
  "$OUT/combined60g" \
  "$TMP" \
  --min-seq-id 0.30 \
  -c 0.70 \
  --cov-mode 0 \
  --cluster-mode 0 \
  -s 7.5 \
  --threads "$THREADS"

python "$HERE/scripts/finalize_combined_family_census.py" \
  --map "$OUT/GLOBAL_EXACT_AA_MAP.tsv" \
  --cluster "$OUT/combined60g_cluster.tsv" \
  --out "$OUT" \
  --mmseqs "$MMSEQS" \
  --mmseqs-version "$MMSEQS_VERSION"

sha256sum \
  "$OUT/COMBINED_GLOBAL_EXACT_UNIQUE.faa" \
  "$OUT/GLOBAL_EXACT_AA_MAP.tsv" \
  "$OUT/combined60g_cluster.tsv" \
  "$OUT/COMBINED_FAMILY_SUMMARY.tsv" \
  "$OUT/STAGE04A_60G_SATURATION_CENSUS.json" \
  > "$OUT/SHA256SUMS.txt"

cat "$OUT/STAGE04A_60G_SATURATION_CENSUS.json"
echo "PASS_STAGE05F3V2_STAGE04A_60G_SATURATION_CENSUS_COMPLETE"
