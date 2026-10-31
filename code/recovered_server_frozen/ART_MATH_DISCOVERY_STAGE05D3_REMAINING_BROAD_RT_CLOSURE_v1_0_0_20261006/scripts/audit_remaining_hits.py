#!/usr/bin/env python3
import csv,re,sys,json
from pathlib import Path
from Bio import SeqIO

D1=Path(sys.argv[1]); ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs'; OUT.mkdir(exist_ok=True)
targets=list(csv.DictReader(open(ROOT/'config/FROZEN_REMAINING_BROAD_RT_TARGETS.tsv'),delimiter='\t'))

# ORF registry
reg=list(csv.DictReader(open(D1/'config/NOVELTY_FIRST_24020_ORF_REGISTRY.tsv'),delimiter='\t'))
bykey={r.get('record_key') or r.get('rt_record_key') or r.get('orf_record_key') or r.get('orf_id'):r for r in reg}
# fallback: first column is commonly record key
if not all(t['rt_record_key'] in bykey for t in targets):
    bykey={next(iter(r.values())):r for r in reg}

# protein fanout maps ORF to unique AA target id/hash when available
fan=list(csv.DictReader(open(D1/'config/NOVELTY_FIRST_PROTEIN_FANOUT.tsv'),delimiter='\t'))
fan_by_orf={}
for r in fan:
    vals=list(r.values())
    joined='\t'.join(str(x) for x in vals)
    for t in targets:
        if t['rt_record_key'] in joined:
            fan_by_orf[t['rt_record_key']]=r

# formal HMM tables: parse non-comment whitespace records
# hmmscan tblout: target model, accession, query seq, accession, E, score...
def tbl_rows(p):
    out=[]
    for line in open(p):
        if not line.strip() or line.startswith('#'): continue
        f=line.split()
        if len(f)>=6: out.append(f)
    return out
T=tbl_rows(D1/'outputs/d1a_frozen/RVT_ALL_E1e4.tblout')
D=tbl_rows(D1/'outputs/d1a_frozen/RVT_ALL_E1e4.domtblout')

# unique proteins FASTA indexed by id
faa={r.id:str(r.seq) for r in SeqIO.parse(D1/'inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa','fasta')}

def infer_aa_id(orf):
    r=fan_by_orf.get(orf,{})
    for v in r.values():
        if isinstance(v,str) and v.startswith('AA_'): return v
    # scan fanout text for AA_ token
    m=re.search(r'AA_[A-Za-z0-9]+','\t'.join(map(str,r.values())))
    return m.group(0) if m else ''

def motif_scan(seq):
    hits=[]
    for pat,name in [(r'Y.DD','YxDD'),(r'[YF].DD','[YF]xDD'),(r'D.{0,3}D','D-x0-3-D')]:
        for m in re.finditer(pat,seq): hits.append(f'{name}:{m.start()+1}-{m.end()}:{m.group(0)}')
    return ';'.join(hits) if hits else 'NONE'

rows=[]
for t in targets:
    orf=t['rt_record_key']; aa_id=infer_aa_id(orf); seq=faa.get(aa_id,'')
    tbl=[f for f in T if aa_id and aa_id in f]
    dom=[f for f in D if aa_id and aa_id in f]
    fullE=''; model=t['myrt_model']
    if tbl:
        # best row by first parsable E after first 4 fields
        best=tbl[0]; model=best[0]
        try: fullE=str(float(best[4]))
        except: fullE=''
    dom_cE=dom_iE=ali=''
    if dom:
        best=dom[0]
        # HMMER domtblout standard columns: c-E 11, i-E 12, hmmfrom 15 hmmto16 alifrom17 alito18 (0-based 11/12/15..18 after split)
        try:
            dom_cE=best[11]; dom_iE=best[12]; ali=f'{best[17]}-{best[18]}'
        except: pass
    rr=bykey.get(orf,{})
    rows.append({
      **t,
      'aa_id':aa_id,
      'observed_aa_len':len(seq) if seq else '',
      'aa_sequence_recovered':bool(seq),
      'formal_model':model,
      'formal_full_E':fullE,
      'formal_domain_cE':dom_cE,
      'formal_domain_iE':dom_iE,
      'formal_alignment_aa':ali,
      'motif_scan':motif_scan(seq) if seq else 'NOT_RECOVERED',
      'orf_registry_source':'NOVELTY_FIRST_24020_ORF_REGISTRY.tsv',
    })

p=OUT/'REMAINING_BROAD_RT_HIT_AUDIT.tsv'
with open(p,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
print('\n'.join('\t'.join(str(r[k]) for k in rows[0]) for r in []))
print(f'[OK] wrote {p} rows={len(rows)}')
