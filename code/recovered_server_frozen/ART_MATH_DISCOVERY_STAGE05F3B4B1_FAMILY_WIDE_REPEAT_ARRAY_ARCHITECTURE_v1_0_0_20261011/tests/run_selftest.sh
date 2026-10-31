#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
python - "$HERE" <<'PY'
import importlib.util,sys
from pathlib import Path
here=Path(sys.argv[1]); p=here/'scripts/repeat_array_scan.py'
spec=importlib.util.spec_from_file_location('b',p); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
motif='ACGTGCTAACGT'
sp=['A'*150,'C'*160,'G'*140]
seq='TTT'+motif+sp[0]+motif+sp[1]+motif+sp[2]+motif+'GGG'
a=m.scan_region(seq,'UPSTREAM',-len(seq))
assert a, 'no array'
best=sorted(a,key=lambda x:(-x['copy_n'],-x['motif_len'],x['spacer_cv']))[0]
assert best['copy_n']>=4,(best,a[:3])
assert best['motif'].find(motif[:8])>=0 or motif.find(best['motif'][:8])>=0,(best['motif'],motif)
assert best['broad_art_like_geometry'] is True,best
# downstream same sequence must never be strict published ART geometry
b=m.scan_region(seq,'DOWNSTREAM',1000)
assert all(not x['strict_published_art_geometry'] for x in b)
print('[PASS] F3B4B1 synthetic repeat-array / orientation / ART-geometry selftest')
PY
