#!/usr/bin/env bash
set -euo pipefail
R4=${1:?R4 package dir containing references/af_swissprot and references/uniprot_sprot.fasta.gz}
THREADS=${2:-32}
mkdir -p work/p1_cif_flat raw outputs logs tmp/full3d_search
rm -f work/p1_cif_flat/*
for d in boltz_out/predictions/*; do
  [ -d "$d" ] || continue
  q=$(basename "$d")
  c=$(find "$d" -maxdepth 1 -name '*_model_0.cif' | head -1)
  [ -n "$c" ] && ln -sfn "$(readlink -f "$c")" "work/p1_cif_flat/${q}.cif"
done
n=$(find work/p1_cif_flat -type l | wc -l); [ "$n" -eq 43 ] || { echo "[FAIL] expected 43 flat CIFs, got $n"; exit 3; }
foldseek easy-search work/p1_cif_flat "$R4/references/af_swissprot" raw/P1_FULL3D_AFSP.m8 tmp/full3d_search -e 0.01 --max-seqs 50 --threads "$THREADS" --exact-tmscore 1 --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,tmscore 2>&1 | tee logs/foldseek_full3d.log
python scripts/summarize_r5.py --freeze frozen --boltz boltz_out --m8 raw/P1_FULL3D_AFSP.m8 --sprot "$R4/references/uniprot_sprot.fasta.gz" --out outputs
cat outputs/GATE_STAGE05D4R5.json
