#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
command -v foldseek >/dev/null || { echo '[ERROR] foldseek not found'; exit 2; }
mkdir -p references tmp/ref_downloads
if [ ! -e references/prostt5 ]; then
  echo '[FETCH] Foldseek ProstT5 weights'
  foldseek databases ProstT5 references/prostt5 tmp/ref_downloads
fi
if [ ! -e references/af_swissprot ]; then
  echo '[FETCH] Foldseek AlphaFold/Swiss-Prot database'
  foldseek databases Alphafold/Swiss-Prot references/af_swissprot tmp/ref_downloads
fi
if [ ! -s references/uniprot_sprot.fasta.gz ]; then
  echo '[FETCH] UniProt Swiss-Prot FASTA for target descriptions'
  curl -fL --retry 4 --retry-delay 3 \
    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
    -o references/uniprot_sprot.fasta.gz
fi
sha256sum references/uniprot_sprot.fasta.gz > references/uniprot_sprot.fasta.gz.sha256
printf '[PASS] references prepared\n'
