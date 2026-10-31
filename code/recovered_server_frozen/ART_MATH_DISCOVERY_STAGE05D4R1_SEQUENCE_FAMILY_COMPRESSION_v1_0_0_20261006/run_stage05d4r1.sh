#!/usr/bin/env bash
set -euo pipefail
D4="${1:?usage: ./run_stage05d4r1.sh /path/to/Stage05D4}"
mkdir -p work/nr90 work/fam50 outputs logs
FA="$D4/outputs/STAGE05D4_N1N2_NONFORMAL_BROAD_RT_PROTEINS.faa"
DB=work/queryDB
rm -f "$DB"* 2>/dev/null || true
mmseqs createdb "$FA" "$DB"
# Primary redundancy collapse: near-identical proteins only.
mmseqs cluster "$DB" work/nr90/clu work/nr90/tmp --min-seq-id 0.90 -c 0.80 --cov-mode 0 --threads "${THREADS:-16}"
mmseqs createtsv "$DB" "$DB" work/nr90/clu work/nr90/nr90_cluster.tsv
# Descriptive family map: broader sequence grouping; never a filter.
mmseqs cluster "$DB" work/fam50/clu work/fam50/tmp --min-seq-id 0.50 -c 0.70 --cov-mode 0 --threads "${THREADS:-16}"
mmseqs createtsv "$DB" "$DB" work/fam50/clu work/fam50/fam50_cluster.tsv
python scripts/build_cluster_tables.py --d4 "$D4" --work work --out outputs | tee logs/run_stage05d4r1.log
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
