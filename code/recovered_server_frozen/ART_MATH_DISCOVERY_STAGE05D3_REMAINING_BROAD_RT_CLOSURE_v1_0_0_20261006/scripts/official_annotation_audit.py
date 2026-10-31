#!/usr/bin/env python3
import csv,sys
from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs'; GB=ROOT/'inputs/official_pc'
targets=list(csv.DictReader(open(ROOT/'config/FROZEN_REMAINING_BROAD_RT_TARGETS.tsv'),delimiter='\t'))
D1=Path(sys.argv[1])
# recover Prodigal proteins directly from C1-like D1 registry isn't sequence-bearing, so use fanout -> unique fasta
fan=list(csv.DictReader(open(D1/'config/NOVELTY_FIRST_PROTEIN_FANOUT.tsv'),delimiter='\t'))
faa={r.id:str(r.seq) for r in SeqIO.parse(D1/'inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa','fasta')}

def aa_for_orf(orf):
    for r in fan:
        s='\t'.join(map(str,r.values()))
        if orf in s:
            for v in r.values():
                if isinstance(v,str) and v.startswith('AA_') and v in faa: return v,faa[v]
    return '', ''

rows=[]
for t in targets:
    acc=t['accession']; gb=GB/f'{acc}.gb'
    base={'accession':acc,'rt_record_key':t['rt_record_key']}
    if not gb.exists():
        rows.append({**base,'status':'OFFICIAL_INPUT_MISSING'}); continue
    rec=SeqIO.read(gb,'genbank')
    aa_id,pseq=aa_for_orf(t['rt_record_key'])
    ps,pe=int(t['rt_start0']),int(t['rt_end0'])
    candidates=[]
    for feat in rec.features:
        if feat.type!='CDS': continue
        s=int(feat.location.start); e=int(feat.location.end); strand='+' if feat.location.strand==1 else '-' if feat.location.strand==-1 else '.'
        ov=max(0,min(pe,e)-max(ps,s)); dist=0 if ov else min(abs(ps-e),abs(s-pe))
        trans=(feat.qualifiers.get('translation') or [''])[0]
        exact=(bool(pseq) and trans==pseq)
        candidates.append((exact,ov,-dist,s,e,strand,feat,trans))
    if not candidates:
        rows.append({**base,'status':'NO_OFFICIAL_CDS'}); continue
    candidates.sort(reverse=True,key=lambda x:(x[0],x[1],x[2]))
    exact,ov,negdist,s,e,strand,feat,trans=candidates[0]
    product=(feat.qualifiers.get('product') or [''])[0]
    pid=(feat.qualifiers.get('protein_id') or [''])[0]
    frac=''
    if pseq and trans:
        al=PairwiseAligner(); al.mode='global'; al.match_score=1; al.mismatch_score=0; al.open_gap_score=0; al.extend_gap_score=0
        frac=float(al.score(pseq,trans))/max(len(pseq),len(trans))
    rows.append({**base,'status':'PASS','aa_id':aa_id,'official_start0':s,'official_end0':e,'official_strand':strand,'official_overlap_bp':ov,'official_distance_bp':-negdist,'official_product':product,'official_protein_id':pid,'pred_vs_official_globalxx_fraction':frac,'exact_translation_match':exact})

cols=[]
for r in rows:
    for k in r:
        if k not in cols: cols.append(k)
with open(OUT/'REMAINING_OFFICIAL_ANNOTATION_AUDIT.tsv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
for r in rows: print('\t'.join(str(r.get(k,'')) for k in cols))
