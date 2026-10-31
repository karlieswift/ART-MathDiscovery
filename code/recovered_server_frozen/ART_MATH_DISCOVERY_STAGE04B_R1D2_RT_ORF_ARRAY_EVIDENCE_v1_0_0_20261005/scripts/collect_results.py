#!/usr/bin/env python3
from pathlib import Path
import tarfile,time,hashlib,sys
from common import *
ROOT=Path(__file__).resolve().parents[1]
if not (ROOT/'outputs').exists():raise SystemExit('No execution outputs to collect')
# Small reproducible inputs are already in the delivery package; do not pack
# 174 upstream archives or any UG27 feature arrays. Include raw ORF and HMM logs.
files=[]
for dirname in ('outputs','logs','refs','config','docs','scripts','tests'):
    for p in (ROOT/dirname).rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts:files.append(p)
for p in (ROOT/'work/prodigal').rglob('*'):
    if p.is_file():files.append(p)
for p in (ROOT/'work/tool_smoke').rglob('*'):
    if p.is_file():files.append(p)
files.extend([ROOT/'README.md',ROOT/'PACKAGE_MANIFEST.tsv',ROOT/'PACKAGE_MANIFEST.sha256',ROOT/'environment.yml'])
files.extend(ROOT.glob('*.sh'))
files=sorted(set(p for p in files if p.is_file()))
manifest=ROOT/'RETURN_MANIFEST.tsv';output_manifest(ROOT,files,manifest);files.append(manifest)
name=ROOT.name+'_RESULTS_'+time.strftime('%Y%m%d_%H%M%S')+'.tar.gz';dest=ROOT/name;tmp=dest.with_suffix('.gz.partial')
if dest.exists():raise SystemExit('Result archive already exists')
with tarfile.open(tmp,'w:gz') as t:
    for p in files:t.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
with tarfile.open(tmp) as t:
    if len(t.getmembers())!=len(files):raise RuntimeError('Archive member count mismatch')
tmp.replace(dest);atomic_text(str(dest)+'.sha256',sha_file(dest)+'  '+dest.name+'\n')
print(dest);print(sha_file(dest));print('bytes',dest.stat().st_size)
