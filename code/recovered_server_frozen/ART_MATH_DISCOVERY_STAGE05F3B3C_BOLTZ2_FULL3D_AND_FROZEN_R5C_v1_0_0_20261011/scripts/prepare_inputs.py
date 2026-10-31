#!/usr/bin/env python3
import csv, hashlib, json, sys
from pathlib import Path

qfa=Path(sys.argv[1]); audit=Path(sys.argv[2]); outdir=Path(sys.argv[3]); manifest=Path(sys.argv[4])
outdir.mkdir(parents=True, exist_ok=True)

def fasta(path):
    h=None; seq=[]
    with path.open() as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: yield h, ''.join(seq)
                h=line[1:].split()[0]; seq=[]
            else: seq.append(line)
        if h is not None: yield h, ''.join(seq)

seqs={h:s for h,s in fasta(qfa)}
if len(seqs)!=32: raise SystemExit(f'FAIL handoff FASTA n={len(seqs)} expected=32')
if len(set(seqs))!=32: raise SystemExit('FAIL duplicate FASTA IDs')

with audit.open() as f:
    rows=list(csv.DictReader(f, delimiter='\t'))
if len(rows)!=32: raise SystemExit(f'FAIL audit rows={len(rows)} expected=32')

cols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','yaml_path','sequence_sha256']
out=[]
seen=set()
for r in rows:
    cid=r['candidate_id']
    if cid in seen: raise SystemExit(f'FAIL duplicate candidate_id {cid}')
    seen.add(cid)
    if cid not in seqs: raise SystemExit(f'FAIL candidate missing from FASTA: {cid}')
    seq=seqs[cid]
    if int(r['aa_length']) != len(seq): raise SystemExit(f'FAIL length mismatch {cid}: {len(seq)} != {r["aa_length"]}')
    y=outdir/f'{cid}.yaml'
    # Boltz v1 YAML schema; MSA intentionally omitted because primary uses --use_msa_server.
    y.write_text('version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: '+seq+'\n')
    out.append({
        'shortlist_rank':int(r['shortlist_rank']), 'family_rep_id':r['family_rep_id'],
        'candidate_id':cid, 'accession':r['accession'], 'aa_length':len(seq),
        'yaml_path':str(y), 'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest()
    })

out.sort(key=lambda x:x['shortlist_rank'])
with manifest.open('w', newline='') as f:
    w=csv.DictWriter(f, delimiter='\t', fieldnames=cols); w.writeheader(); w.writerows(out)
print('query_n =',len(out))
print('max_aa_length =',max(x['aa_length'] for x in out))
print('PASS_F3B3C_INPUT_PREPARATION')
