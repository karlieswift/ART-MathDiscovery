#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
python - "$S/scripts/reconstruct_survivor_loci.py" <<'PY'
import importlib.util, pathlib, sys
modp=pathlib.Path(sys.argv[1]).resolve()
spec=importlib.util.spec_from_file_location('m',modp);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
pep='MKTWQ'
seq='CCC'+'ATGAAAACTTGGCAA'+'TAA'+'GGG'
h=m.exact_peptide_matches(seq,pep)
assert h==[('+',0,3,18)],h
seq2='AAAA'+m.revcomp('ATGAAAACTTGGCAA')+'CCCC'
h2=m.exact_peptide_matches(seq2,pep)
assert len(h2)==1 and h2[0][0]=='-' and h2[0][2]==4 and h2[0][3]==19,h2
assert m.translate_frame('ATGAAATAA',0)=='MK*'
print('[PASS] F3B4A exact-placement selftest')
PY
