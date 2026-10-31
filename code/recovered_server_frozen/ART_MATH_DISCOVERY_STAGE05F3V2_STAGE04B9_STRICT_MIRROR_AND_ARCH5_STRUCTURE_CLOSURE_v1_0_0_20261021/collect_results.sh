#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd);STAMP=$(date +%Y%m%d_%H%M%S);OUT="$HERE/../ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B9_RESULTS_${STAMP}.tar.gz"
cd "$HERE"
FILES=(README.md config/CONTRACT.json docs/CLAIM_BOUNDARY.md docs/SUPPLEMENT_ROUTE_TABLE_TEMPLATE.md)
for f in outputs/*.json outputs/*.tsv outputs/*.txt logs/*.log; do [[ -e "$f" ]] && FILES+=("$f"); done
tar -czf "$OUT" "${FILES[@]}"
sha256sum "$OUT" > "$OUT.sha256"
echo "$OUT";cat "$OUT.sha256"
