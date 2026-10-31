#!/usr/bin/env bash
set -euo pipefail
P=${1:?Usage: bash run_n1_external_screen.sh PROJECT_ROOT [--metadata /path/to/metadata.tsv] [--gbff-dir /path/to/gbff]}
shift
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$P/exports/N1_EXTERNAL_SCREEN_$(date +%Y%m%d_%H%M%S)_$$"
mkdir -p "$OUT"
python "$ROOT/scripts/curate_public_sources.py" --kit "$ROOT" --out "$OUT"
while (($#)); do
 case "$1" in
 --metadata)
   test $# -ge 2;python "$ROOT/scripts/screen_new_phage_metadata.py" --metadata "$2" --exposure "$ROOT/data/EXPOSURE_919_AND_HISTORICAL_ACCESSIONS.tsv" --out "$OUT";shift 2;;
 --gbff-dir)
   test $# -ge 2;python "$ROOT/scripts/screen_genbank_rt_features.py" --gbff-dir "$2" --exposure "$ROOT/data/EXPOSURE_919_AND_HISTORICAL_ACCESSIONS.tsv" --out "$OUT";shift 2;;
 *) echo "Unknown argument: $1" >&2;exit 2;;
 esac
done
ARCHIVE="$OUT.tar.gz"
tar -czf "$ARCHIVE" -C "$(dirname "$OUT")" "$(basename "$OUT")"
sha256sum "$ARCHIVE" > "$ARCHIVE.sha256"
sha256sum -c "$ARCHIVE.sha256"
printf 'N1_STATUS=%s\nSEND_ARCHIVE=%s\nSEND_SHA256=%s\n' "$OUT/N1_PUBLIC_CURATION_GATE.json" "$ARCHIVE" "$ARCHIVE.sha256"
