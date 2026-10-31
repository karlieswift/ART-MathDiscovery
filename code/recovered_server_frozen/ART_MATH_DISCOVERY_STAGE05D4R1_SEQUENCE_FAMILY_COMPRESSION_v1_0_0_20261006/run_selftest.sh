#!/usr/bin/env bash
set -euo pipefail
D4="${1:?usage: ./run_selftest.sh /path/to/Stage05D4}"
python tests/selftest.py "$D4"
command -v mmseqs >/dev/null || { echo '[FAIL] mmseqs not found. Install once: conda install -y -c conda-forge -c bioconda mmseqs2'; exit 2; }
mmseqs version
