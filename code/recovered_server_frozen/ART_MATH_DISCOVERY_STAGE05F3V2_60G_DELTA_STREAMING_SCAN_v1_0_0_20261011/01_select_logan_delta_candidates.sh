#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd); P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
K=${1:-5000}
mkdir -p "$S"/{inputs,logs,frozen}
[[ -s "$S/inputs/OLD_F3_ACCESSIONS.txt" ]] || python "$S/scripts/discover_old_accessions.py" "$P" "$S/inputs/OLD_F3_ACCESSIONS.txt" | tee "$S/logs/discover_old_accessions.log"
N=$(grep -cve '^$' "$S/inputs/OLD_F3_ACCESSIONS.txt" || true)
if (( N < 131 )); then
  echo "[FATAL] OLD_F3_ACCESSIONS.txt has $N accessions; require >=131 before freezing a truly delta-only universe." >&2
  echo "Add the missing old F3 source accessions to $S/inputs/OLD_F3_ACCESSIONS.txt and rerun." >&2
  exit 2
fi
sha256sum "$S/inputs/OLD_F3_ACCESSIONS.txt" > "$S/frozen/OLD_F3_ACCESSIONS.sha256"
echo "[STREAM] Logan metadata -> deterministic source-blind candidate pool (K=$K/stratum)"
aws s3 cp s3://logan-pub/stats/sra_taxid.csv.zst - --no-sign-request \
 | zstd -dc \
 | python "$S/scripts/select_logan_metadata.py" --exclude "$S/inputs/OLD_F3_ACCESSIONS.txt" --k "$K" --out "$S/inputs/F3V2_DELTA_ACCESSION_CANDIDATES.tsv" \
 2> >(tee "$S/logs/select_logan_candidates.stderr" >&2)
sha256sum "$S/inputs/F3V2_DELTA_ACCESSION_CANDIDATES.tsv" > "$S/frozen/F3V2_DELTA_CANDIDATE_POOL.sha256"
echo PASS_F3V2_DELTA_SOURCE_BLIND_CANDIDATE_POOL_FROZEN
