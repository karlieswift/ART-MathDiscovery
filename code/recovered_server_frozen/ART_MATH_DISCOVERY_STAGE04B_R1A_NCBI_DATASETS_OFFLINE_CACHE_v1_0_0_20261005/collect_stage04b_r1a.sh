#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
cat outputs/STAGE04B_R1A_GATE.json
wc -l outputs/STAGE04B_R1A_DATASETS_CACHE_MANIFEST.tsv
sha256sum outputs/STAGE04B_R1A_DATASETS_CACHE_MANIFEST.tsv > outputs/STAGE04B_R1A_DATASETS_CACHE_MANIFEST.tsv.sha256
cat outputs/STAGE04B_R1A_DATASETS_CACHE_MANIFEST.tsv.sha256
