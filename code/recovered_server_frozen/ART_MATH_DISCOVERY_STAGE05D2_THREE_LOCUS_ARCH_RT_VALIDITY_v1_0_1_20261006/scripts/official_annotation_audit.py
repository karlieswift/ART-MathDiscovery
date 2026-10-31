#!/usr/bin/env python3
import sys, json
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
D=Path(sys.argv[1]); F=pd.read_csv(sys.argv[2],sep='\t'); gbdir=Path(sys.argv[3]); out=Path(sys.argv[4])
fan=pd.read_csv(D/'config/NOVELTY_FIRST_PROTEIN_FANOUT.tsv',sep='\t')
seqs={r.id:str(r.seq) for r in SeqIO.parse(D/'inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa','fasta')}
rows=[]
for _,f in F.iterrows():
    acc=f.accession; target=f.broad_rt_target
    pred=fan[(fan.sequence_id==acc)&(fan.target_id==target)].iloc[0]
    g=gbdir/f'{acc}.gb'
    if not g.exists():
        rows.append(dict(accession=acc,candidate_key=f.candidate_key,status='OFFICIAL_FETCH_FAILED'))
        continue
    try:
        rec=SeqIO.read(g,'genbank')
    except Exception as e:
        rows.append(dict(accession=acc,candidate_key=f.candidate_key,status='OFFICIAL_PARSE_FAILED',error_type=type(e).__name__,error_message=str(e)[:240]))
        continue
    best=None
    for feat in rec.features:
        if feat.type!='CDS': continue
        s=int(feat.location.start); e=int(feat.location.end); strand=feat.location.strand
        ov=max(0,min(e,int(pred.end0))-max(s,int(pred.start0)))
        dist=0 if ov else min(abs(s-int(pred.end0)),abs(int(pred.start0)-e))
        tr=(feat.qualifiers.get('translation') or [''])[0]
        if tr:
            aligner=PairwiseAligner()
            aligner.mode='global'
            aligner.match_score=1.0
            aligner.mismatch_score=0.0
            aligner.open_gap_score=0.0
            aligner.extend_gap_score=0.0
            a=aligner.score(seqs[target],tr)
            ident=a/max(len(seqs[target]),len(tr))
        else: ident=0.0
        key=(ov,ident,-dist)
        if best is None or key>best[0]: best=(key,s,e,strand,feat,ident,ov,dist)
    if best is None:
        rows.append(dict(accession=acc,candidate_key=f.candidate_key,status='NO_OFFICIAL_CDS'))
        continue
    _,s,e,strand,feat,ident,ov,dist=best
    rows.append(dict(accession=acc,candidate_key=f.candidate_key,status='PASS',predicted_record_key=pred.record_key,pred_start0=int(pred.start0),pred_end0=int(pred.end0),pred_strand=pred.strand,
                     official_start0=s,official_end0=e,official_strand='+' if strand==1 else '-' if strand==-1 else '.',official_overlap_bp=ov,official_distance_bp=dist,
                     official_product=(feat.qualifiers.get('product') or [''])[0],official_protein_id=(feat.qualifiers.get('protein_id') or [''])[0],pred_vs_official_globalxx_fraction=ident))
pd.DataFrame(rows).to_csv(out/'OFFICIAL_ANNOTATION_AUDIT.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
