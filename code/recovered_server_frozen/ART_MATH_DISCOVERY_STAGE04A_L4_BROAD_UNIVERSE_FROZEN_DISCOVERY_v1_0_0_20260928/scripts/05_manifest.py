#!/usr/bin/env python3
from pathlib import Path
import argparse, hashlib
from common import sha256_file

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); a=ap.parse_args(); root=Path(a.root).resolve(); out=root/'STAGE04A_PACKAGE_MANIFEST.tsv'
    exclude={out.resolve()}
    rows=[]
    for p in sorted(root.rglob('*')):
        if not p.is_file() or p.resolve() in exclude: continue
        if 'runtime_engine' in p.parts or 'work' in p.parts: continue
        rows.append((str(p.relative_to(root)),p.stat().st_size,sha256_file(p)))
    with open(out,'w') as f:
        f.write('path\tsize_bytes\tsha256\n')
        for r in rows: f.write(f'{r[0]}\t{r[1]}\t{r[2]}\n')
    print(f'[OK] manifest files={len(rows)} sha256={sha256_file(out)}')
if __name__=='__main__': main()
