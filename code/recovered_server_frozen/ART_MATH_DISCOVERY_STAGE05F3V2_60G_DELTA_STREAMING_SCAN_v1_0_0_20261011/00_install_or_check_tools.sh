#!/usr/bin/env bash
set -euo pipefail
ENV=${CONDA_ENV_NAME:-art_novelty}
missing=()
for x in aws zstd hmmscan hmmpress mmseqs python sha256sum; do command -v "$x" >/dev/null 2>&1 || missing+=("$x"); done
if (( ${#missing[@]} )); then
  echo "[MISSING] ${missing[*]}"
  echo "Install once with:"
  echo "conda install -n $ENV -y -c conda-forge -c bioconda awscli zstd hmmer mmseqs2"
  exit 2
fi
echo '[PASS] aws zstd HMMER MMseqs2 Python available'
