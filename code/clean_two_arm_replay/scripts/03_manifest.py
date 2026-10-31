#!/usr/bin/env python3
import argparse,hashlib
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);a=ap.parse_args();root=Path(a.root).resolve();rows=[]
for p in sorted(root.rglob('*')):
    if p.is_file() and 'logs/' not in str(p.relative_to(root)) and p.name!='PACKAGE_MANIFEST.sha256':
        h=hashlib.sha256(p.read_bytes()).hexdigest();rows.append((h,str(p.relative_to(root))))
(root/'PACKAGE_MANIFEST.sha256').write_text(''.join(f'{h}  {n}\n' for h,n in rows));print(f'[OK] manifest files={len(rows)}')
