#!/usr/bin/env python3
import argparse,csv,gzip,hashlib,json,os,re

def fasta_records(path):
    op=gzip.open if path.endswith('.gz') else open
    name=None; seq=[]
    with op(path,'rt',errors='replace') as f:
        for ln in f:
            if ln.startswith('>'):
                if name is not None: yield name,''.join(seq)
                name=ln[1:].rstrip(); seq=[]
            else: seq.append(ln.strip())
        if name is not None: yield name,''.join(seq)

def parse_header(h):
    tok=h.split()[0]; acc=tok; entry=''
    if '|' in tok:
        p=tok.split('|')
        if len(p)>=3: acc,entry=p[1],p[2]
        elif len(p)>=2: acc=p[1]
    desc=h[len(tok):].strip()
    return acc,entry,desc

def parse_lookup(path):
    byacc={}; rx=re.compile(r'AF-([A-Za-z0-9]+)-F\d+-model',re.I)
    with open(path,errors='replace') as f:
        for ln in f:
            if not ln.strip(): continue
            p=ln.rstrip('\n').split('\t')
            if len(p)<2: continue
            key,name=p[0],p[1]
            m=rx.search(name)
            if m: byacc.setdefault(m.group(1),[]).append((key,name))
    return byacc

BAD_TERMS=[
'reverse transcriptase','rna-directed dna polymerase','telomerase','retron',
'group ii intron','group 2 intron','intron-encoded protein','maturase',
'retrotranspos','retroelement','retrovirus','retroviral','ribonuclease h-pol',
]

def safe_negative(desc):
    d=desc.lower()
    return not any(t in d for t in BAD_TERMS)

def hkey(seed,pos,acc):
    return hashlib.sha256((seed+'|'+pos+'|'+acc).encode()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--panel',required=True); ap.add_argument('--sprot',required=True); ap.add_argument('--lookup',required=True)
    ap.add_argument('--outdir',required=True); ap.add_argument('--per-rt',type=int,default=3); ap.add_argument('--seed',default='ART_STAGE05D4R5C_20261007')
    a=ap.parse_args(); os.makedirs(a.outdir,exist_ok=True)
    panel=list(csv.DictReader(open(a.panel),delimiter='\t'))
    panel_acc={r['accession'] for r in panel}
    clean=[r for r in panel if r.get('annotation_class') and r.get('has_afsp_structure')=='true' and r.get('afsp_key')]
    lookup=parse_lookup(a.lookup)
    meta={}
    for hdr,seq in fasta_records(a.sprot):
        acc,entry,desc=parse_header(hdr); meta[acc]=(entry,desc,len(seq))
    pos=[]
    for r in clean:
        acc=r['accession']; entry,desc,L=meta[acc]
        pos.append({'accession':acc,'entry':entry,'description':desc,'aa_length':L,'annotation_class':r['annotation_class'],'afsp_key':r['afsp_key'],'afsp_target':r['afsp_target']})
    # Candidate negative pool: AF structure available, outside panel, no RT-associated annotation text.
    negpool=[]
    for acc,(entry,desc,L) in meta.items():
        if acc in panel_acc or acc not in lookup or not safe_negative(desc): continue
        # choose the first deterministic AF mapping (normally one)
        key,name=sorted(lookup[acc],key=lambda z:(int(z[0]) if z[0].isdigit() else 10**18,z[0]))[0]
        negpool.append((acc,entry,desc,L,key,name))
    used=set(); neg=[]; tols=[10,25,50,100,200]
    for p in sorted(pos,key=lambda x:x['accession']):
        chosen=[]
        for tol in tols:
            cands=[x for x in negpool if x[0] not in used and abs(x[3]-p['aa_length'])<=tol]
            cands.sort(key=lambda x:(abs(x[3]-p['aa_length']),hkey(a.seed,p['accession'],x[0]),x[0]))
            for x in cands:
                if len(chosen)>=a.per_rt: break
                used.add(x[0]); chosen.append(x)
            if len(chosen)>=a.per_rt: break
        if len(chosen)<a.per_rt:
            raise SystemExit(f'Could not find {a.per_rt} unique negative controls for {p["accession"]}; got {len(chosen)}')
        for x in chosen:
            acc,entry,desc,L,key,name=x
            neg.append({'matched_rt_accession':p['accession'],'matched_rt_length':p['aa_length'],'accession':acc,'entry':entry,'description':desc,'aa_length':L,'length_delta':L-p['aa_length'],'afsp_key':key,'afsp_target':name})
    def write(path,rows,fields):
        with open(path,'w',newline='') as f:
            w=csv.DictWriter(f,delimiter='\t',fieldnames=fields);w.writeheader();w.writerows(rows)
    write(os.path.join(a.outdir,'STAGE05D4R5C_CLEAN_RT_REFERENCES.tsv'),pos,['accession','entry','description','aa_length','annotation_class','afsp_key','afsp_target'])
    write(os.path.join(a.outdir,'STAGE05D4R5C_MATCHED_NONRT_CONTROLS.tsv'),neg,['matched_rt_accession','matched_rt_length','accession','entry','description','aa_length','length_delta','afsp_key','afsp_target'])
    with open(os.path.join(a.outdir,'STAGE05D4R5C_CLEAN_RT_KEYS.txt'),'w') as f:
        for x in sorted({r['afsp_key'] for r in pos}, key=lambda s:(int(s) if s.isdigit() else 10**18,s)): f.write(x+'\n')
    with open(os.path.join(a.outdir,'STAGE05D4R5C_NONRT_KEYS.txt'),'w') as f:
        for x in sorted({r['afsp_key'] for r in neg}, key=lambda s:(int(s) if s.isdigit() else 10**18,s)): f.write(x+'\n')
    S={'stage':'Stage05D4-R5C-freeze','n_clean_rt_structures':len(pos),'n_matched_nonrt_controls':len(neg),'negative_controls_per_rt':a.per_rt,'seed':a.seed,'reference_selection_uses_candidate_results':False,'gate':'PASS_CLEAN_RT_AND_MATCHED_NONRT_CALIBRATION_UNIVERSE_FROZEN' if len(pos)>=30 and len(neg)==len(pos)*a.per_rt else 'FAIL_CALIBRATION_FREEZE','gate_pass':len(pos)>=30 and len(neg)==len(pos)*a.per_rt}
    json.dump(S,open(os.path.join(a.outdir,'STAGE05D4R5C_FREEZE.json'),'w'),indent=2)
    print(json.dumps(S,indent=2))
if __name__=='__main__': main()
