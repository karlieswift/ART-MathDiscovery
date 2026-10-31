#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
from common import *
KEYS=['OR836606','LC680885','MN091626','MZ779063','CP001363','CP026202','AP006725']
TOK=['truth','array','repeat','locus','rt_start','rt_end','array_start','array_end','core_start','core_end']
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('project_root'); a=ap.parse_args(); root=Path(a.project_root); here=Path(__file__).resolve().parents[1]
    rows=[]; hints=('STAGE03B','STAGE04C1','STAGE05C2','ART_PANEL','ART95','REVIEW_OR770614','TRUTH','UNBLIND')
    for p in root.rglob('*'):
        if not p.is_file() or p.suffix.lower() not in {'.tsv','.csv','.json','.txt'}:continue
        if here in p.parents:continue
        up=str(p).upper()
        if not any(h in up for h in hints):continue
        try:
            if p.stat().st_size>20_000_000:continue
            head=p.read_text(errors='replace')[:120000]
        except Exception:continue
        hu=head.upper()
        if not any(k in up or k in hu for k in KEYS):continue
        low=head.lower(); score=sum(t in low for t in TOK)
        if score:rows.append({'score':score,'accessions':';'.join(k for k in KEYS if k in up or k in hu),'path':str(p),'sha256':sha256_file(p)})
    rows.sort(key=lambda r:(-r['score'],r['path'])); out=here/'outputs/TRUTH_SOURCE_CANDIDATES.tsv'
    if rows:write_tsv(out,rows,['score','accessions','path','sha256'])
    print(f'wrote {out} rows={len(rows)}')
if __name__=='__main__':main()
