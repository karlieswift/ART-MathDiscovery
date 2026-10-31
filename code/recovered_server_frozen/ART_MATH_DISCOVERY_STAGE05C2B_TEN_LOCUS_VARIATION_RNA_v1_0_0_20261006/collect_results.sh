#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p outputs work
exec 9>work/.run.lock
flock -n 9 || { echo "[HOLD] Stop or complete the active run before collecting results."; exit 2; }
TS=$(date +%Y%m%d_%H%M%S)
NAME="ART_STAGE05C2B_TEN_LOCUS_RESULTS_${TS}.tar.gz"
python - <<'PY'
from pathlib import Path
import hashlib,csv
base=Path('.');rows=[]
for d in ['outputs','logs','work/fold','work/author']:
 if not Path(d).exists():continue
 for p in sorted(Path(d).rglob('*')):
  if not p.is_file() or p.is_symlink() or p.name=='RESULT_MANIFEST.tsv' or p.name=='technical_report_text_local.txt' or p.name.endswith('.part'):continue
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(1<<20),b''):h.update(b)
  rows.append([str(p),p.stat().st_size,h.hexdigest()])
with open('outputs/RESULT_MANIFEST.tsv','w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['path','size_bytes','sha256']);w.writerows(rows)
PY
ITEMS=(outputs inputs precomputed config scripts tests docs README.md PACKAGE_MANIFEST.tsv PACKAGE_MANIFEST.sha256 run_selftest.sh run_stage05c2b.sh collect_results.sh fetch_author_sources.sh reproduce_precomputed.sh)
for X in logs work/fold work/author; do if [ -e "$X" ]; then ITEMS+=("$X"); fi; done
# Full source PDF is included when acquired so the returned record can be checked.
# Extracted full article text and unfinished transfers are not included.
tar --exclude='technical_report_text_local.txt' --exclude='*.part' --exclude='__pycache__' -czf "${NAME}.part" "${ITEMS[@]}"
gzip -t "${NAME}.part"
mv "${NAME}.part" "$NAME"
sha256sum "$NAME" > "${NAME}.sha256"
sha256sum -c "${NAME}.sha256"
ls -lh "$NAME" "${NAME}.sha256"
