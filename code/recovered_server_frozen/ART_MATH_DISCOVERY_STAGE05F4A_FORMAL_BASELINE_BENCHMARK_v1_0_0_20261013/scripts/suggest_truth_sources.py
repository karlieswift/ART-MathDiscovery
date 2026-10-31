#!/usr/bin/env python3
from __future__ import annotations
import argparse,re
from pathlib import Path
from common import *
KEYS=['OR836606','LC680885','MN091626','MZ779063','NC_016856','CP026202']
TOK=['truth','array','repeat','locus','rt_start','rt_end','array_start','array_end','core_start','core_end']
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('project_root'); a=ap.parse_args(); root=Path(a.project_root)
    rows=[]
    hints=('STAGE03B','STAGE04C1','STAGE05C2','ART_PANEL','ART95','REVIEW_OR770614','TRUTH','UNBLIND')
    for p in root.rglob('*'):
        if not p.is_file() or p.suffix.lower() not in {'.tsv','.csv','.json','.txt'}: continue
        up=str(p).upper()
        if not any(h in up for h in hints): continue
        try:
            if p.stat().st_size > 20_000_000: continue
        except Exception:
            continue
        u=str(p).upper()
        if not any(k in u for k in KEYS):
            try: head=p.read_text(errors='replace')[:120000]
            except: continue
            hu=head.upper()
            if not any(k in hu for k in KEYS): continue
        else:
            try: head=p.read_text(errors='replace')[:120000]
            except: continue
        low=head.lower(); score=sum(t in low for t in TOK)
        if score:
            present=';'.join(k for k in KEYS if k in head.upper())
            rows.append({'score':score,'accessions':present,'path':str(p),'sha256':sha256_file(p)})
    rows.sort(key=lambda r:(-r['score'],r['path']))
    out=Path(__file__).resolve().parents[1]/'outputs/TRUTH_SOURCE_CANDIDATES.tsv'
    if rows: write_tsv(out,rows,['score','accessions','path','sha256'])
    print(f'wrote {out} rows={len(rows)}')
if __name__=='__main__': main()
