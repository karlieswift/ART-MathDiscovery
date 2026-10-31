#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
python - "$HERE" <<'PY'
import importlib.util,sys
from pathlib import Path
here=Path(sys.argv[1])
def load(name,p):
 s=importlib.util.spec_from_file_location(name,p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
r=load('r',here/'scripts/vendor_repeat_array_scan.py')
p=load('p',here/'scripts/vendor_family_wide_context_rescue.py')
q=load('q',here/'scripts/partner_orf_scan.py')
# Frozen repeat detector: 3 exact 8-mer copies with 150-nt spacers => strict ART geometry.
m='ACGTGCAA'; spacer='C'*150
seq='G'*20 + m + spacer + m + spacer + m + 'T'*20
arr=r.scan_region(seq,'UPSTREAM',-len(seq))
assert any(a['strict_published_art_geometry'] for a in arr), arr
# Exact peptide placement.
dna='TAA'+'ATG'*40+'TAA'
pep='M'*40
h=p.exact_hits(dna,pep)
assert len(h)==1 and h[0][0]=='+',h
# Partner ORF caller: complete same-strand ORF >=30 aa.
partner='ATG'+'GCT'*59+'TAA'
ors=q.call_forward_orfs(partner,p.CODON,'SAME_STRAND')
assert len(ors)==1 and ors[0]['aa_length']==60,ors
print('[PASS] F3B4C frozen detector / exact-placement / partner-ORF synthetic selftest')
PY
