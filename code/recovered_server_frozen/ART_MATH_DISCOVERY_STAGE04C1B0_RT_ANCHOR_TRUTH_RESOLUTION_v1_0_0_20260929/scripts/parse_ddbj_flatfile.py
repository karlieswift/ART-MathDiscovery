#!/usr/bin/env python3
import re, json, argparse
from pathlib import Path

def parse_loc(s):
    s=s.strip(); strand='+'
    if s.startswith('complement(') and s.endswith(')'):
        strand='-'; s=s[len('complement('):-1]
    # join is allowed; use outer span
    nums=[int(x) for x in re.findall(r'\d+',s)]
    if not nums: return None
    return min(nums),max(nums),strand

def parse_features(path):
    feats=[]; cur=None; lastq=None
    for raw in Path(path).read_text(errors='replace').splitlines():
        if raw.startswith('     ') and len(raw)>=21 and raw[5:21].strip() and not raw[21:].lstrip().startswith('/'):
            key=raw[5:21].strip(); loc=raw[21:].strip()
            if key and re.match(r'^[A-Za-z_]+$',key):
                if cur: feats.append(cur)
                p=parse_loc(loc)
                cur={'type':key,'location':loc,'start1':p[0] if p else None,'end1':p[1] if p else None,'strand':p[2] if p else None,'qual':{}}
                lastq=None
                continue
        if cur and raw.startswith('                     /'):
            x=raw.strip()[1:]
            if '=' in x:
                k,v=x.split('=',1); v=v.strip().strip('"'); cur['qual'].setdefault(k,[]).append(v); lastq=k
            else:
                cur['qual'].setdefault(x,[]).append(''); lastq=x
        elif cur and lastq and raw.startswith('                     ') and raw.strip() and not raw.strip().startswith('/'):
            # continuation of quoted qualifier
            cur['qual'][lastq][-1]+=raw.strip().strip('"')
    if cur: feats.append(cur)
    return feats

def q(f,k): return ';'.join(f['qual'].get(k,[]))
def aa_len(f):
    tr=q(f,'translation').replace(' ','')
    return len(tr) if tr else None

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--flatfile',required=True); ap.add_argument('--out',required=True)
    a=ap.parse_args(); fs=parse_features(a.flatfile)
    rows=[]
    for f in fs:
        if f['type'] not in {'CDS','gene','ncRNA','misc_RNA','regulatory'}: continue
        rows.append({**{k:f[k] for k in ('type','start1','end1','strand','location')},'locus_tag':q(f,'locus_tag'),'gene':q(f,'gene'),'product':q(f,'product'),'protein_id':q(f,'protein_id'),'aa_len':aa_len(f)})
    Path(a.out).write_text(json.dumps(rows,indent=2))
if __name__=='__main__': main()
