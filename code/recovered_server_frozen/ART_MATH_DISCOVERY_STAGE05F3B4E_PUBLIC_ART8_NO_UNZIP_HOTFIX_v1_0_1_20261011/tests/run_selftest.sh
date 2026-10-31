#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
bash -n "$S/run_stage05f3b4e_public_art8_r1.sh"
python - <<'PY'
import zipfile,tempfile,os
p=tempfile.mktemp('.zip')
with zipfile.ZipFile(p,'w') as z:z.writestr('ncbi_dataset/data/genomic.fna','>MW218148.1 x\nACGT\n')
with zipfile.ZipFile(p) as z: assert 'ncbi_dataset/data/genomic.fna' in z.namelist()
os.unlink(p)
print('[PASS] Python zipfile available')
PY
echo '[PASS] F3B4E public ART8 R1 syntax/extraction selftest'
