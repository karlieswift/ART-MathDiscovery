#!/usr/bin/env bash
set -euo pipefail
D4="${1:?usage: ./apply_hotfix_and_finish.sh /path/to/Stage05D4 [Stage05D4R1_v1_0_0_dir]}"
R1="${2:-$PWD}"
FIXDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$R1"
for f in work/nr90/nr90_cluster.tsv work/fam50/fam50_cluster.tsv; do
  test -s "$f" || { echo "[FAIL] missing completed MMseqs cluster table: $f" >&2; exit 2; }
done
cp -f scripts/build_cluster_tables.py scripts/build_cluster_tables.py.v1_0_0.bak
cp -f "$FIXDIR/scripts/build_cluster_tables.py" scripts/build_cluster_tables.py
rm -f outputs/STAGE05D4R1_* outputs/GATE_STAGE05D4R1.json 2>/dev/null || true
python scripts/build_cluster_tables.py --d4 "$D4" --work work --out outputs | tee logs/run_stage05d4r1_hotfix_v1_0_1.log
python - <<'PY'
from pathlib import Path
import hashlib
out=Path('outputs'); rows=[]
for p in sorted(out.iterdir()):
    if p.is_file() and p.name!='STAGE05D4R1_OUTPUT_MANIFEST.tsv':
        h=hashlib.sha256(p.read_bytes()).hexdigest(); rows.append((p.as_posix(),p.stat().st_size,h))
with open(out/'STAGE05D4R1_OUTPUT_MANIFEST.tsv','w') as f:
    f.write('file\tbytes\tsha256\n')
    for r in rows:f.write('\t'.join(map(str,r))+'\n')
PY
echo '[PASS] Stage05D4-R1 v1.0.1 namespace hotfix applied; existing MMseqs clustering reused unchanged.'
