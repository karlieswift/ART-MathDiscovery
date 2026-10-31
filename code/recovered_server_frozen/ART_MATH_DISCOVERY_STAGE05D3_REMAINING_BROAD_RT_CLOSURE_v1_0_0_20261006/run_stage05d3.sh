#!/usr/bin/env bash
set -euo pipefail
D1=${1:?usage: ./run_stage05d3.sh /path/to/Stage05D1}
mkdir -p outputs
for acc in MT708547.1 NC_031049.1 PV251870.1; do
  f="inputs/official_pc/${acc}.gb"
  if [ ! -s "$f" ]; then
    echo "[FAIL] missing $f; retrieve on Windows PC using windows_download_official_genbank.ps1" >&2
    exit 2
  fi
  grep -q '^LOCUS' "$f" || { echo "[FAIL] invalid GenBank: $f" >&2; exit 2; }
  grep -q '^ORIGIN' "$f" || { echo "[FAIL] no ORIGIN: $f" >&2; exit 2; }
done
python scripts/audit_remaining_hits.py "$D1"
python scripts/official_annotation_audit.py "$D1"
sha256sum inputs/official_pc/*.gb > outputs/OFFICIAL_PC_GENBANK_SHA256.tsv
python scripts/write_gate.py
python scripts/write_report.py
cat outputs/GATE_STAGE05D3.json
