#!/usr/bin/env python3
import csv, re, sys, json, os
from collections import defaultdict

m8, priority, fasta, outdir, sprot_fa = sys.argv[1:6]
os.makedirs(outdir, exist_ok=True)

# query metadata from R3
meta={}
with open(priority, newline='') as f:
    r=csv.DictReader(f, delimiter='\t')
    for row in r:
        if row.get('structure_priority_tier')=='TIER_A_COMPACT_STRUCTURE_READY':
            meta[row['query']] = row

# fasta IDs, for exact universe check
fids=[]
with open(fasta) as f:
    for line in f:
        if line.startswith('>'):
            fids.append(line[1:].split()[0])
if set(fids)!=set(meta):
    raise SystemExit(f"Tier-A FASTA/priority mismatch: fasta={len(fids)} meta={len(meta)}")

# SwissProt descriptions. Accept sp|ACC|ENTRY, ACC, or AF-ACC-F1-like target ids.
desc={}
if sprot_fa!='NA' and os.path.exists(sprot_fa):
    op = open
    if sprot_fa.endswith('.gz'):
        import gzip
        op = gzip.open
    with op(sprot_fa, 'rt') as f:
        for line in f:
            if not line.startswith('>'): continue
            h=line[1:].rstrip()
            tok=h.split()[0]
            acc=tok
            if '|' in tok:
                p=tok.split('|')
                if len(p)>=2: acc=p[1]
            description=h[len(tok):].strip()
            desc[acc]=description

rt_re=re.compile(r'(reverse transcriptase|rna-directed dna polymerase|telomerase reverse transcriptase|retr(?:on|oelement).*reverse transcriptase|group ii intron.*(?:maturase|reverse transcriptase))', re.I)

def normalize_target(t):
    x=t
    if x.startswith('AF-') and '-F' in x:
        x=x[3:].split('-F',1)[0]
    if '|' in x:
        p=x.split('|')
        if len(p)>=2: x=p[1]
    return x

hits=defaultdict(list)
if os.path.exists(m8) and os.path.getsize(m8)>0:
    with open(m8, newline='') as f:
        r=csv.DictReader(f, delimiter='\t', fieldnames=['query','target','evalue','bits','qcov','tcov','alnlen','qlen','tlen','fident'])
        for row in r:
            if row['query'] not in meta: continue
            acc=normalize_target(row['target'])
            d=desc.get(acc,'')
            row['target_accession']=acc
            row['target_description']=d
            row['rt_text_hit']=bool(rt_re.search(d)) if d else False
            hits[row['query']].append(row)

allhits=[]; summary=[]
for q in fids:
    hs=hits.get(q,[])
    # preserve Foldseek order, but robustly choose numeric best evalue
    def ev(h):
        try:return float(h['evalue'])
        except:return float('inf')
    best=min(hs,key=ev) if hs else None
    rths=[h for h in hs if h['rt_text_hit']]
    bestrt=min(rths,key=ev) if rths else None
    m=meta[q]
    label='NO_AFSP_STRUCTURAL_HIT_NOT_NEGATIVE'
    if hs: label='CURATED_STRUCTURAL_NON_RT_TRIAGE'
    if rths: label='STRUCTURAL_RT_TEXT_HIT_REQUIRES_VALIDATION'
    summary.append({
        'query':q,
        'aa_length':m.get('aa_length',''),
        'candidate_key_n':m.get('candidate_key_n_y',m.get('candidate_key_n_x','')),
        'candidate_keys':m.get('arch_candidate_keys',m.get('candidate_keys','')),
        'strata':m.get('strata_y',m.get('strata_x','')),
        'min_orf_to_candidate_distance_bp':m.get('min_orf_to_candidate_distance_bp',''),
        'n_structural_hits':len(hs),
        'best_target':best['target'] if best else '',
        'best_target_accession':best['target_accession'] if best else '',
        'best_evalue':best['evalue'] if best else '',
        'best_bits':best['bits'] if best else '',
        'best_qcov':best['qcov'] if best else '',
        'best_tcov':best['tcov'] if best else '',
        'best_description':best['target_description'] if best else '',
        'n_rt_text_hits':len(rths),
        'best_rt_target':bestrt['target'] if bestrt else '',
        'best_rt_accession':bestrt['target_accession'] if bestrt else '',
        'best_rt_evalue':bestrt['evalue'] if bestrt else '',
        'best_rt_qcov':bestrt['qcov'] if bestrt else '',
        'best_rt_tcov':bestrt['tcov'] if bestrt else '',
        'best_rt_description':bestrt['target_description'] if bestrt else '',
        'triage_label':label
    })
    for h in hs:
        z={'query':q, **h}
        allhits.append(z)

sum_path=os.path.join(outdir,'STAGE05D4R4_TIER_A_STRUCTURAL_TRIAGE.tsv')
fields=list(summary[0].keys())
with open(sum_path,'w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=fields); w.writeheader(); w.writerows(summary)

hit_path=os.path.join(outdir,'STAGE05D4R4_ALL_STRUCTURAL_HITS.tsv')
if allhits:
    hfields=list(allhits[0].keys())
    with open(hit_path,'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=hfields); w.writeheader(); w.writerows(allhits)
else:
    open(hit_path,'w').write('query\ttarget\n')

counts=defaultdict(int)
for x in summary: counts[x['triage_label']]+=1
S={
 'stage':'Stage05D4-R4',
 'n_tier_A_queries':len(summary),
 'n_queries_with_any_structural_hit':sum(x['n_structural_hits']>0 for x in summary),
 'n_queries_with_rt_text_structural_hit':sum(x['n_rt_text_hits']>0 for x in summary),
 'triage_counts':dict(counts),
 'swissprot_description_mapping_available':bool(desc)
}
with open(os.path.join(outdir,'STAGE05D4R4_SUMMARY.json'),'w') as f: json.dump(S,f,indent=2)
print(json.dumps(S,indent=2))
