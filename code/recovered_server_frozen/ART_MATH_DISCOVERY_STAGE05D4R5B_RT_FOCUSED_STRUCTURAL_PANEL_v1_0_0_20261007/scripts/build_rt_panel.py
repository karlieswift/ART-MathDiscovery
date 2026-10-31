#!/usr/bin/env python3
import argparse,csv,gzip,json,os,re,hashlib,collections

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def parse_fasta_headers(path):
    op=gzip.open if path.endswith('.gz') else open
    rows=[]
    with op(path,'rt') as f:
        for line in f:
            if not line.startswith('>'): continue
            h=line[1:].rstrip()
            tok=h.split()[0]
            acc=tok
            entry=''
            if '|' in tok:
                p=tok.split('|')
                if len(p)>=3:
                    acc=p[1]; entry=p[2]
                elif len(p)>=2: acc=p[1]
            desc=h[len(tok):].strip()
            rows.append((acc,entry,desc))
    return rows

def annotation_class(desc):
    d=desc.lower()
    if 'telomerase reverse transcriptase' in d: return 'TELOMERASE_RT'
    if 'group ii intron' in d and ('reverse transcriptase' in d or 'maturase' in d): return 'GROUP_II_INTRON_RT'
    if 'retron' in d and ('reverse transcriptase' in d or 'rna-directed dna polymerase' in d): return 'RETRON_RT'
    if ('diversity-generating' in d or 'diversity generating' in d or 'dgr' in d) and ('reverse transcriptase' in d or 'rna-directed dna polymerase' in d): return 'DGR_RT'
    if ('crispr' in d or 'cas1' in d) and ('reverse transcriptase' in d or 'rna-directed dna polymerase' in d): return 'CRISPR_ASSOCIATED_RT'
    if 'reverse transcriptase' in d or 'rna-directed dna polymerase' in d: return 'GENERIC_RT_TEXT'
    return ''

def parse_hmmtbl(path, topn):
    # hmmsearch tblout: target name/accession, query name/accession, full E-value, score, ...
    per=collections.defaultdict(list)
    if not os.path.exists(path): return per
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith('#'): continue
            x=ln.split()
            if len(x)<6: continue
            target=x[0]; model=x[2]
            try:e=float(x[4]); score=float(x[5])
            except:continue
            acc=target
            if '|' in target:
                p=target.split('|')
                if len(p)>=2: acc=p[1]
            per[model].append((e,-score,acc,target,score))
    selected=collections.defaultdict(list)
    for model,arr in per.items():
        seen=set()
        for e,ns,acc,target,score in sorted(arr):
            if acc in seen: continue
            seen.add(acc); selected[model].append((acc,e,score))
            if len(selected[model])>=topn: break
    return selected

def parse_lookup(path):
    byacc={}
    bykey={}
    afre=re.compile(r'AF-([A-Za-z0-9]+)-F\d+-model',re.I)
    with open(path,errors='replace') as f:
        for ln in f:
            if not ln.strip(): continue
            p=ln.rstrip('\n').split('\t')
            if len(p)<2: continue
            key=p[0]; name=p[1]
            m=afre.search(name)
            if not m: continue
            acc=m.group(1)
            byacc.setdefault(acc,[]).append((key,name))
            bykey[key]=(acc,name)
    return byacc,bykey

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--sprot',required=True)
    ap.add_argument('--hmmtbl',required=True)
    ap.add_argument('--lookup',required=True)
    ap.add_argument('--top-per-profile',type=int,default=10)
    ap.add_argument('--outdir',required=True)
    a=ap.parse_args(); os.makedirs(a.outdir,exist_ok=True)
    headers=parse_fasta_headers(a.sprot)
    ann={}
    meta={}
    for acc,entry,desc in headers:
        c=annotation_class(desc)
        meta[acc]=(entry,desc)
        if c: ann[acc]=c
    hmm=parse_hmmtbl(a.hmmtbl,a.top_per_profile)
    hacc=collections.defaultdict(list); he=collections.defaultdict(list)
    for model,arr in hmm.items():
        for acc,e,score in arr:
            hacc[acc].append(model); he[acc].append((e,model,score))
    union=sorted(set(ann)|set(hacc))
    byacc,_=parse_lookup(a.lookup)
    rows=[]; keys=[]
    for acc in union:
        entry,desc=meta.get(acc,('',''))
        models=sorted(set(hacc.get(acc,[])))
        best=min(he.get(acc,[(float('inf'),'','')]), key=lambda z:z[0])
        maps=byacc.get(acc,[])
        if not maps:
            rows.append(dict(accession=acc,entry=entry,description=desc,annotation_class=ann.get(acc,''),source_arms=';'.join([x for x in ['ANNOTATION' if acc in ann else '', 'RVT_ALL_HMM_TOP' if acc in hacc else ''] if x]),hmm_models=';'.join(models),best_hmm_evalue='' if best[0]==float('inf') else best[0],best_hmm_model=best[1],afsp_key='',afsp_target='',has_afsp_structure='false'))
            continue
        # AlphaFold/Swiss-Prot should normally have one entry; retain all mappings deterministically.
        for key,name in sorted(maps,key=lambda z:(int(z[0]) if z[0].isdigit() else 10**18,z[0])):
            keys.append(key)
            rows.append(dict(accession=acc,entry=entry,description=desc,annotation_class=ann.get(acc,''),source_arms=';'.join([x for x in ['ANNOTATION' if acc in ann else '', 'RVT_ALL_HMM_TOP' if acc in hacc else ''] if x]),hmm_models=';'.join(models),best_hmm_evalue='' if best[0]==float('inf') else best[0],best_hmm_model=best[1],afsp_key=key,afsp_target=name,has_afsp_structure='true'))
    tsv=os.path.join(a.outdir,'STAGE05D4R5B_RT_REFERENCE_PANEL.tsv')
    fields=['accession','entry','description','annotation_class','source_arms','hmm_models','best_hmm_evalue','best_hmm_model','afsp_key','afsp_target','has_afsp_structure']
    with open(tsv,'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields); w.writeheader(); w.writerows(rows)
    keys=sorted(set(keys),key=lambda x:(int(x) if x.isdigit() else 10**18,x))
    with open(os.path.join(a.outdir,'STAGE05D4R5B_RT_REFERENCE_KEYS.txt'),'w') as f:
        for k in keys:f.write(k+'\n')
    prof=[]
    for model in sorted(hmm):
        accs=[x[0] for x in hmm[model]]
        mapped=sum(1 for x in accs if x in byacc)
        prof.append({'hmm_model':model,'selected_swissprot_n':len(accs),'afsp_mapped_n':mapped,'accessions':';'.join(accs)})
    with open(os.path.join(a.outdir,'STAGE05D4R5B_PROFILE_COVERAGE.tsv'),'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=['hmm_model','selected_swissprot_n','afsp_mapped_n','accessions']);w.writeheader();w.writerows(prof)
    S={
      'stage':'Stage05D4-R5B-reference-freeze',
      'swissprot_entries_total':len(headers),
      'annotation_rt_accessions':len(ann),
      'hmm_profiles_with_hits':len(hmm),
      'hmm_top_union_accessions':len(hacc),
      'union_reference_accessions':len(union),
      'afsp_structural_reference_entries':len(keys),
      'reference_panel_selected_from_candidate_results':False,
      'top_per_hmm_profile':a.top_per_profile,
      'gate':'PASS_RT_FOCUSED_REFERENCE_PANEL_FROZEN' if len(keys)>=20 else 'FAIL_RT_REFERENCE_PANEL_TOO_SMALL',
      'gate_pass':len(keys)>=20
    }
    with open(os.path.join(a.outdir,'STAGE05D4R5B_REFERENCE_FREEZE.json'),'w') as f:json.dump(S,f,indent=2)
    print(json.dumps(S,indent=2))
if __name__=='__main__':main()
