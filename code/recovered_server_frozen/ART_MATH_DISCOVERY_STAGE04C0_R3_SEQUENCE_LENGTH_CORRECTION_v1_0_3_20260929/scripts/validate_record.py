#!/usr/bin/env python3
import argparse, hashlib, json
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
    if not seq: raise ValueError('empty FASTA sequence')
    bad=set(seq.upper())-set('ACGTNRYKMSWBDHVU-')
    if bad: raise ValueError('non-IUPAC sequence characters: '+''.join(sorted(bad))[:80])
    return len(seq), lines[0]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--control-id',required=True); ap.add_argument('--registry-accession',required=True); ap.add_argument('--fetch-accession',required=True)
    ap.add_argument('--relation',required=True); ap.add_argument('--expected-length-bp',type=int,required=True)
    ap.add_argument('--fasta',required=True); ap.add_argument('--flatfile',required=True); ap.add_argument('--out',required=True)
    a=ap.parse_args(); fa=Path(a.fasta); gb=Path(a.flatfile)
    for p in (fa,gb):
        if not p.is_file() or p.stat().st_size==0: raise SystemExit(f'[FAIL] missing/empty {p}')
        head=p.read_text(errors='replace')[:1000].lower()
        if '<!doctype html' in head or '<html' in head or 'www error blocked diagnostic' in head:
            raise SystemExit(f'[FAIL] HTML/block page detected: {p}')
        if p.stat().st_size < 1000: raise SystemExit(f'[FAIL] suspiciously tiny record: {p}')
    n,header=fasta_len(fa)
    if n != a.expected_length_bp:
        raise SystemExit(f'[FAIL] unexpected sequence length {n} != frozen expected {a.expected_length_bp}')
    fetch0=a.fetch_accession.split('.')[0]
    if fetch0 not in header: raise SystemExit(f'[FAIL] accession {fetch0} absent FASTA header: {header}')
    gbtext=gb.read_text(errors='replace')
    if not gbtext.startswith('LOCUS'): raise SystemExit('[FAIL] flatfile does not start with LOCUS')
    if fetch0 not in gbtext[:10000]: raise SystemExit(f'[FAIL] accession {fetch0} absent flatfile header')
    out={
      'gate':'PASS_STAGE04C0_R3_VALID_RECORD','gate_pass':True,'control_id':a.control_id,
      'registry_accession':a.registry_accession,'accession':a.fetch_accession,'fetch_accession':a.fetch_accession,'relation':a.relation,
      'length_bp':n,'expected_length_bp':a.expected_length_bp,'fasta_header':header,
      'fasta_sha256':sha256(fa),'flatfile_sha256':sha256(gb),'fasta_bytes':fa.stat().st_size,'flatfile_bytes':gb.stat().st_size
    }
    Path(a.out).write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
