#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$(pwd)}"
mkdir -p "$ROOT/inputs" "$ROOT/logs"
XLSX="$ROOT/inputs/Minerva_SupplementaryTables_preprint.xlsx"
PDF="$ROOT/inputs/Minerva_manuscript_joined.pdf"
XLSX_URL='https://zenodo.org/records/22903368/files/Minerva_SupplementaryTables_preprint.xlsx?download=1'
PDF_URL='https://zenodo.org/records/22903368/files/Minerva_manuscript_joined.pdf?download=1'
fetch() {
  local url="$1" out="$2"
  if command -v wget >/dev/null 2>&1; then
    wget -c --tries=3 --timeout=30 "$url" -O "$out"
  elif command -v curl >/dev/null 2>&1; then
    curl -L --retry 3 --connect-timeout 30 -o "$out" "$url"
  else
    echo "[ERROR] neither wget nor curl found" >&2; return 2
  fi
}
fetch "$XLSX_URL" "$XLSX"
# PDF is useful for cross-checking figure/method wording but is not required for the resolver gate.
fetch "$PDF_URL" "$PDF" || true
sha256sum "$XLSX" > "$ROOT/inputs/Minerva_SupplementaryTables_preprint.xlsx.sha256"
if [[ -s "$PDF" ]]; then sha256sum "$PDF" > "$ROOT/inputs/Minerva_manuscript_joined.pdf.sha256"; fi
ls -lh "$ROOT/inputs" | tee "$ROOT/logs/download_sources.log"
