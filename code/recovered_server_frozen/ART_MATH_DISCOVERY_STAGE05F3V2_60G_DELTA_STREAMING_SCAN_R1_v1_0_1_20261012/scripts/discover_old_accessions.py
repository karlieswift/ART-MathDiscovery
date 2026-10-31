#!/usr/bin/env python3
import os,re,sys,csv
from pathlib import Path
P=Path(sys.argv[1]); out=Path(sys.argv[2])
pat=re.compile(r'\b([SED]RR\d+)\b')
acc=set(); sources=[]
# First: compact project tables/manifests likely to contain the old accessions.
for p in P.rglob('*'):
    try:
        if not p.is_file(): continue
        name=p.name.lower(); sp=str(p)
        if p.stat().st_size>100_000_000: continue
        if not any(k in name for k in ['manifest','instance_map','source','accession','gate']): continue
        if 'f3' not in sp.lower(): continue
        if p.suffix.lower() not in ['.tsv','.csv','.txt','.json','.jsonl']: continue
        txt=p.read_text(errors='ignore')
        got=set(pat.findall(txt))
        if got:
            acc |= got; sources.append((str(p),len(got)))
    except Exception:
        pass
# If we already have >=131, write. We only need to ensure all old source accessions are excluded;
# over-exclusion is scientifically harmless for the delta universe.
if len(acc)<131:
    # Search a few large F3 FASTA candidates and inspect headers only.
    cands=[]
    for p in P.rglob('*'):
        try:
            if not p.is_file(): continue
            if p.suffix.lower() not in ['.fa','.fna','.fasta','.fas']: continue
            if p.stat().st_size < 500_000_000: continue
            if 'f3' not in str(p).lower(): continue
            cands.append(p)
        except Exception: pass
    cands=sorted(cands,key=lambda p:p.stat().st_size,reverse=True)[:12]
    for p in cands:
        got=set()
        try:
            with open(p,errors='ignore') as f:
                for line in f:
                    if line.startswith('>'):
                        got.update(pat.findall(line))
            if got:
                acc |= got; sources.append((str(p),len(got)))
            if len(acc)>=131: break
        except Exception: pass
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(''.join(x+'\n' for x in sorted(acc)))
print('old_exclusion_n =',len(acc))
for p,n in sources[-20:]: print('SOURCE',n,p)
if len(acc)<98:
    raise SystemExit('FATAL: could not recover even the 98 formal-source accessions; supply inputs/OLD_F3_ACCESSIONS.txt manually')
if len(acc)<131:
    print('WARNING: recovered <131 old source accessions; before final source freeze, add the remaining old F3 source accessions to inputs/OLD_F3_ACCESSIONS.txt',file=sys.stderr)
