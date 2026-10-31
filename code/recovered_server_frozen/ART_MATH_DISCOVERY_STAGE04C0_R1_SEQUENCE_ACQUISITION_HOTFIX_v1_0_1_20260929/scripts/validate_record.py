#!/usr/bin/env python3
import argparse, hashlib, json, re
from pathlib import Path

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def fasta_len(p):
    lines=Path(p).read_text(errors='replace').splitlines()
    if not lines or not lines[0].startswith('>'):
        raise ValueError('FASTA does not start with >')
    seq=''.join(x.strip() for x in lines[1:] if x.strip())
    if not seq:
        raise ValueError('empty FASTA sequence')
    bad=set(seq.upper())-set('ACGTNRYKMSWBDHVU-')
    if bad:
        raise ValueError('non-IUPAC sequence characters: '+''.join(sorted(bad))[:80])
    return len(seq), lines[0]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--control-id',required=True); ap.add_argument('--accession',required=True)
    ap.add_argument('--fasta',required=True); ap.add_argument('--flatfile',required=True)
    ap.add_argument('--min-bp',type=int,default=1000000); ap.add_argument('--out',required=True)
    a=ap.parse_args()
    fa=Path(a.fasta); gb=Path(a.flatfile)
    for p in (fa,gb):
        if not p.is_file() or p.stat().st_size==0: raise SystemExit(f'[FAIL] missing/empty {p}')
        head=p.read_text(errors='replace')[:500]
        if '<!DOCTYPE html' in head or '<html' in head.lower() or 'WWW Error Blocked Diagnostic' in head:
            raise SystemExit(f'[FAIL] HTML/block page detected: {p}')
    n,header=fasta_len(fa)
    if n < a.min_bp: raise SystemExit(f'[FAIL] sequence too short {n} < {a.min_bp}')
    gbtext=gb.read_text(errors='replace')
    if not gbtext.startswith('LOCUS'):
        raise SystemExit('[FAIL] flatfile does not start with LOCUS')
    # DDBJ/INSDC flatfile should contain accession token without requiring exact formatting.
    acc0=a.accession.split('.')[0]
    if acc0 not in gbtext[:5000]: raise SystemExit(f'[FAIL] accession {acc0} absent from flatfile header')
    out={
      'gate':'PASS_STAGE04C0_R1_VALID_RECORD', 'gate_pass':True,
      'control_id':a.control_id,'accession':a.accession,'length_bp':n,
      'fasta_header':header,'fasta_sha256':sha256(fa),'flatfile_sha256':sha256(gb),
      'fasta_bytes':fa.stat().st_size,'flatfile_bytes':gb.stat().st_size
    }
    Path(a.out).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
