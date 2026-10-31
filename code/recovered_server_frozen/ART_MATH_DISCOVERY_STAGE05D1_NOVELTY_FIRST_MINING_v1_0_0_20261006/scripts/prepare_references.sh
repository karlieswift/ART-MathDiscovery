#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REF="$ROOT/references/broad_rt"
mkdir -p "$REF"
fetch() {
  local out="$1"; shift
  if [ -s "$out" ]; then return 0; fi
  for u in "$@"; do
    echo "[FETCH] $u" >&2
    if curl -fL --retry 2 --connect-timeout 20 "$u" -o "$out.tmp"; then mv "$out.tmp" "$out"; return 0; fi
    rm -f "$out.tmp"
  done
  return 1
}
fetch "$REF/RVT-All.hmm" \
  "https://omics.informatics.indiana.edu/myRT/share/RVT-All.hmm" \
  "https://omics.sice.indiana.edu/myRT/share/RVT-All.hmm" \
  "https://raw.githubusercontent.com/mgtools/myRT/master/Models/RVT-All.hmm"
# Full cdd+Pfam refinement database from myRT authors. Large, but only primary hits are queried.
fetch "$REF/cdd-pfamA.tar.gz" \
  "https://omics.informatics.indiana.edu/myRT/share/cdd-pfamA.tar.gz" \
  "https://omics.sice.indiana.edu/myRT/Models/cdd-pfamA.tar.gz"
if [ ! -s "$REF/cdd-pfamA.hmm" ]; then
  tar -xzf "$REF/cdd-pfamA.tar.gz" -C "$REF"
  FOUND=$(find "$REF" -type f -name 'cdd-pfamA.hmm' | head -1 || true)
  [ -n "$FOUND" ] || { echo '[FAIL] cdd-pfamA.hmm not found after extraction' >&2; exit 2; }
  if [ "$FOUND" != "$REF/cdd-pfamA.hmm" ]; then cp "$FOUND" "$REF/cdd-pfamA.hmm"; fi
fi
# Freeze bytes before scanning.
sha256sum "$REF/RVT-All.hmm" "$REF/cdd-pfamA.hmm" > "$REF/REFERENCE_SHA256.txt"
grep -c '^NAME ' "$REF/RVT-All.hmm" | awk '{print "myRT_models=" $1}' > "$REF/REFERENCE_COUNTS.txt"
grep -c '^NAME ' "$REF/cdd-pfamA.hmm" | awk '{print "cdd_pfam_models=" $1}' >> "$REF/REFERENCE_COUNTS.txt"
cat "$REF/REFERENCE_SHA256.txt"
cat "$REF/REFERENCE_COUNTS.txt"
