#!/usr/bin/env python3
import argparse, csv, hashlib, json, re
from pathlib import Path

def read_fasta(p):
    d={}; h=None; s=[]
    with open(p) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: d[h]=''.join(s)
                h=line[1:].split()[0]; s=[]
            else: s.append(line)
    if h is not None: d[h]=''.join(s)
    return d

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def yq(s):
    return '"'+s.replace('\\','\\\\').replace('"','\\"')+'"'

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--r3',required=True); ap.add_argument('--r4',required=True); ap.add_argument('--out',required=True)
    a=ap.parse_args(); r3=Path(a.r3); r4=Path(a.r4); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    g3=json.load(open(r3/'outputs/GATE_STAGE05D4R3.json')); g4=json.load(open(r4/'outputs/GATE_STAGE05D4R4.json'))
    assert g3['gate_pass'] and g3['n_tier_A']==350
    assert g4['gate_pass'] and g4['n_tier_A_queries']==350 and g4['n_queries_with_any_structural_hit']==106 and g4['n_queries_with_rt_text_structural_hit']==0
    rows=[]
    with open(r4/'outputs/STAGE05D4R4_TIER_A_STRUCTURAL_TRIAGE.tsv') as f:
        for x in csv.DictReader(f,delimiter='\t'):
            x['aa_length']=int(x['aa_length']); x['candidate_key_n']=int(x['candidate_key_n']); x['min_orf_to_candidate_distance_bp']=int(float(x['min_orf_to_candidate_distance_bp']))
            rows.append(x)
    assert len(rows)==350
    nohit=[x for x in rows if x['triage_label']=='NO_AFSP_STRUCTURAL_HIT_NOT_NEGATIVE']
    assert len(nohit)==244
    p1=[x for x in nohit if 220<=x['aa_length']<=500 and x['min_orf_to_candidate_distance_bp']<=250 and x['candidate_key_n']>=2]
    p2=[x for x in nohit if 220<=x['aa_length']<=500 and x['min_orf_to_candidate_distance_bp']<=250 and x['candidate_key_n']==1]
    p1.sort(key=lambda x:x['query']); p2.sort(key=lambda x:x['query'])
    assert len(p1)==43, len(p1); assert len(p2)==40, len(p2)
    fa=read_fasta(r3/'outputs/STAGE05D4R3_TIER_A_COMPACT_STRUCTURE_READY.faa')
    assert len(fa)==350
    for x in p1+p2: assert x['query'] in fa
    fields=list(rows[0].keys())+['r5_compute_tier']
    with open(out/'STAGE05D4R5_P1_FREEZE.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader()
        for x in p1: w.writerow({**x,'r5_compute_tier':'P1_REPLICATED_COMPACT_FULL3D'})
    with open(out/'STAGE05D4R5_P2_EXTENSION_FREEZE.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader()
        for x in p2: w.writerow({**x,'r5_compute_tier':'P2_SINGLETON_COMPACT_EXTENSION'})
    with open(out/'STAGE05D4R5_P1.faa','w') as f:
        for x in p1: f.write(f">{x['query']}\n{fa[x['query']]}\n")
    yd=out/'p1_yaml'; yd.mkdir(exist_ok=True)
    for x in p1:
        q=x['query']; seq=fa[q]
        with open(yd/f'{q}.yaml','w') as f:
            f.write('version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: '+yq(seq)+'\n')
    summary={'stage':'Stage05D4-R5-freeze','n_r4_tier_A':350,'n_r4_nohit':244,'n_P1':43,'n_P2':40,'p1_rule':'R4 no-hit; 220-500 aa; min ORF-to-candidate distance <=250 bp; candidate_key_n >=2','p2_rule':'same except candidate_key_n ==1','candidate_filtering_performed':False,'claim_boundary':'P1/P2 are compute-priority tiers only; all R4 unresolved and all 602 DNA-first candidates remain retained.'}
    json.dump(summary,open(out/'STAGE05D4R5_FREEZE_SUMMARY.json','w'),indent=2)
    manifest=[]
    for p in sorted(out.rglob('*')):
        if p.is_file(): manifest.append({'path':str(p.relative_to(out)),'sha256':sha(p),'bytes':p.stat().st_size})
    with open(out/'STAGE05D4R5_FREEZE_MANIFEST.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['path','sha256','bytes'],delimiter='\t'); w.writeheader(); w.writerows(manifest)
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
