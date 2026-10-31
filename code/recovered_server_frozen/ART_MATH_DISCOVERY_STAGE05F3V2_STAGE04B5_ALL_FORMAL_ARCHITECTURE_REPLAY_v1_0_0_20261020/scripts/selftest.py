#!/usr/bin/env python3
import importlib.util, tempfile
from pathlib import Path
p=Path(__file__).with_name('stage04b5.py')
s=importlib.util.spec_from_file_location('b5',p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m)
# exact candidate parsing / translation
x=m.parse_cid('SRR1|ctg1|+|f0|0-300|aa100|o1')
assert x['accession']=='SRR1' and x['start0']==0 and x['end0']==300 and x['aa_length']==100
assert m.translate('ATGGCTTAA')=='MA*'
# strict upstream synthetic array: motif 8, spacer 157 => start gap 165
mot='ATGGAGTA'; arr=mot+'A'*157+mot+'C'*157+mot
seq=arr+'G'*50+'ATG'+'GCT'*99
rt0=len(arr)+50; rt1=rt0+300
z=m.repeat_arrays(seq,rt0,rt1)
assert any(a['strict_published_art_geometry'] for a in z), z
# downstream ORF
seq='A'*300+'ATG'+'GCT'*35+'TAA'+'A'*50
orfs=m.find_orfs_downstream(seq,300)
assert orfs and len(orfs[0][3])>=30
print('PASS_STAGE04B5_SELFTEST')
