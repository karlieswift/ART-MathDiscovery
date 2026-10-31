#!/usr/bin/env python3
import argparse,csv,hashlib,json,os,shutil
from pathlib import Path

def fasta(path):
    d={};h=None;s=[]
    with open(path) as f:
        for ln in f:
            ln=ln.strip()
            if not ln: continue
            if ln.startswith('>'):
                if h is not None:d[h]=''.join(s)
                h=ln[1:].split()[0];s=[]
            else:s.append(ln)
        if h is not None:d[h]=''.join(s)
    return d

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def yq(s): return '"'+s.replace('\\','\\\\').replace('"','\\"')+'"'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--r3',required=True);ap.add_argument('--r5',required=True);ap.add_argument('--out',required=True)
    a=ap.parse_args();r3=Path(a.r3);r5=Path(a.r5);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    p2f=r5/'frozen/STAGE05D4R5_P2_EXTENSION_FREEZE.tsv'
    if not p2f.exists(): raise SystemExit(f'missing {p2f}; run R5 freeze first')
    rows=list(csv.DictReader(open(p2f),delimiter='\t'))
    if len(rows)!=40: raise SystemExit(f'expected 40 P2 rows, got {len(rows)}')
    for r in rows:
        if r.get('r5_compute_tier')!='P2_SINGLETON_COMPACT_EXTENSION' or int(r['candidate_key_n'])!=1:
            raise SystemExit(f'P2 freeze invariant failed for {r.get("query")}')
    fa=fasta(r3/'outputs/STAGE05D4R3_TIER_A_COMPACT_STRUCTURE_READY.faa')
    yd=out/'p2_yaml';yd.mkdir(exist_ok=True)
    with open(out/'STAGE05D4R5D_P2.faa','w') as fo:
        for r in sorted(rows,key=lambda x:x['query']):
            q=r['query'];seq=fa.get(q)
            if not seq: raise SystemExit(f'missing P2 sequence {q} in R3 tier-A FASTA')
            if len(seq)!=int(r['aa_length']): raise SystemExit(f'length mismatch {q}: {len(seq)} vs {r["aa_length"]}')
            fo.write(f'>{q}\n{seq}\n')
            with open(yd/f'{q}.yaml','w') as f:
                f.write('version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: '+yq(seq)+'\n')
    shutil.copy2(p2f,out/'STAGE05D4R5D_P2_FREEZE.tsv')
    s={'stage':'Stage05D4-R5D-freeze','n_P2':40,'source':'Stage05D4-R5 P2 extension freeze','candidate_key_n_rule':1,'candidate_filtering_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'gate':'PASS_P2_CALIBRATED_FULL3D_EXTENSION_FROZEN','gate_pass':True}
    json.dump(s,open(out/'STAGE05D4R5D_FREEZE.json','w'),indent=2)
    with open(out/'STAGE05D4R5D_FREEZE_SHA256.tsv','w') as f:
        for p in sorted(out.rglob('*')):
            if p.is_file() and p.name!='STAGE05D4R5D_FREEZE_SHA256.tsv': f.write(f'{sha(p)}  {p.relative_to(out)}\n')
    print(json.dumps(s,indent=2))
if __name__=='__main__':main()
