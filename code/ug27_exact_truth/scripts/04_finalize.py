#!/usr/bin/env python3
from pathlib import Path
import argparse, json
from common import sha256, dump_json

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--outputs',required=True); a=ap.parse_args(); out=Path(a.outputs)
    files=sorted([p for p in out.rglob('*') if p.is_file() and p.name!='STAGE04C2_R4_SHA256.tsv'])
    lines=[]
    for p in files: lines.append(f"{sha256(p)}\t{p.relative_to(out)}")
    (out/'STAGE04C2_R4_SHA256.tsv').write_text('\n'.join(lines)+'\n')
    print(f'[PASS] hashed {len(files)} R4 output files')
if __name__=='__main__': main()
