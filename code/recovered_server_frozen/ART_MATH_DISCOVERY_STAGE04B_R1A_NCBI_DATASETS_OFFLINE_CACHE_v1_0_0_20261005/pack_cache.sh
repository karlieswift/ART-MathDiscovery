#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
python - <<'PY'
import json
j=json.load(open('outputs/STAGE04B_R1A_GATE.json'))
assert j['gate_pass'] is True, j
assert j['n_cache_rows']==174, j
PY
PKG="STAGE04B_R1A_NCBI_DATASETS_CACHE_174.tar.gz"
tar -czf "$PKG" work/zips work/unpacked outputs/STAGE04B_R1A_DATASETS_CACHE_MANIFEST.tsv outputs/STAGE04B_R1A_DATASETS_CACHE_MANIFEST.tsv.sha256 outputs/STAGE04B_R1A_GATE.json inputs/STAGE04B_R1A_ACCESSIONS_174.txt
sha256sum "$PKG" > "$PKG.sha256"
ls -lh "$PKG" "$PKG.sha256"
cat "$PKG.sha256"
