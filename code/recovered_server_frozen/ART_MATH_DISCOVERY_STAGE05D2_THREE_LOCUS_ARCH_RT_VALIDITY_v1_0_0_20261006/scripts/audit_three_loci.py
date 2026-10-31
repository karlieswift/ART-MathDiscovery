#!/usr/bin/env python3
import sys, re, json, hashlib
from pathlib import Path
import pandas as pd
from Bio import SeqIO
D=Path(sys.argv[1]); frozen=Path(sys.argv[2]); out=Path(sys.argv[3]); out.mkdir(parents=True,exist_ok=True)
F=pd.read_csv(frozen,sep='\t')
C=pd.read_csv(D/'config/NOVELTY_FIRST_999_CANDIDATES.tsv',sep='\t')
L=pd.read_csv(D/'config/NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv',sep='\t')
R=pd.read_csv(D/'config/NOVELTY_FIRST_NONCODING_REPEAT_GROUPS.tsv',sep='\t')
P=pd.read_csv(D/'config/NOVELTY_FIRST_PROTEIN_FANOUT.tsv',sep='\t')
seqs={r.id:str(r.seq) for r in SeqIO.parse(D/'inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa','fasta')}

def parse_dom(path):
    rows=[]
    for line in open(path):
        if not line.strip() or line.startswith('#'): continue
        f=line.split()
        # domtblout: target=model, query=protein sequence
        rows.append(dict(model=f[0], model_len=int(f[2]), target_id=f[3], aa_len=int(f[5]), full_E=float(f[6]), full_score=float(f[7]),
                         i_E=float(f[12]), dom_score=float(f[13]), hmm_from=int(f[15]), hmm_to=int(f[16]), ali_from=int(f[17]), ali_to=int(f[18]), env_from=int(f[19]), env_to=int(f[20]), acc=float(f[21])))
    return pd.DataFrame(rows)
DHM=parse_dom(D/'outputs/d1a_frozen/RVT_ALL_E1e4.domtblout')
rows=[]; reps=[]; motifs=[]
for _,f in F.iterrows():
    ck=f.candidate_key; acc=f.accession; target=f.broad_rt_target
    c=C[C.candidate_key==ck].iloc[0]
    p=P[(P.target_id==target)&(P.sequence_id==acc)]
    if p.empty: raise SystemExit(f'missing target {target} in {acc}')
    p=p.iloc[0]
    l=L[(L.candidate_key==ck)&(L.record_key==p.record_key)]
    if l.empty: raise SystemExit(f'missing link {ck} {p.record_key}')
    l=l.iloc[0]
    rr=R[R.candidate_key==ck].copy()
    # use only the N1/N2-defining prediction-relative noncoding high-complexity groups
    if 'high_complexity' in rr: rr=rr[rr.high_complexity.astype(int)==1]
    if 'full_span_outside_predicted_CDS' in rr: rr=rr[rr.full_span_outside_predicted_CDS.astype(str).isin(['True','1','true'])]
    rep_start=int(rr.start0.min()) if len(rr) else None
    rep_end=int(rr.end0.max()) if len(rr) else None
    if rep_start is None: rep_gap=None
    elif int(p.end0)<=rep_start: rep_gap=rep_start-int(p.end0)
    elif int(p.start0)>=rep_end: rep_gap=int(p.start0)-rep_end
    else: rep_gap=0
    # strand-aware side: genomic lower-coordinate repeat relative to minus ORF = 3' side; plus = 5' side
    if rep_start is None: side='NA'
    elif rep_end <= int(p.start0): side='5prime_upstream' if p.strand=='+' else '3prime_downstream'
    elif rep_start >= int(p.end0): side='3prime_downstream' if p.strand=='+' else '5prime_upstream'
    else: side='overlap'
    d=DHM[(DHM.target_id==target)&(DHM.model==f.myrt_model)]
    if d.empty: raise SystemExit(f'missing dom hit {target} {f.myrt_model}')
    d=d.sort_values(['i_E','dom_score'],ascending=[True,False]).iloc[0]
    model_cov=(d.hmm_to-d.hmm_from+1)/d.model_len
    query_cov=(d.ali_to-d.ali_from+1)/d.aa_len
    aa=seqs[target]
    yxdd=[]
    for m in re.finditer(r'Y.DD',aa): yxdd.append((m.start()+1,m.group()))
    # position relative to significant aligned region
    motif_in_aln=[x for x in yxdd if d.ali_from <= x[0] <= d.ali_to]
    rows.append(dict(candidate_key=ck,accession=acc,region_rank=int(c.region_rank),candidate_start0=int(c.start0),candidate_end0=int(c.end0),stratum=c.novelty_architecture_stratum,
                     rt_record_key=p.record_key,rt_start0=int(p.start0),rt_end0=int(p.end0),rt_strand=p.strand,rt_aa_len=int(p.aa_length),candidate_to_rt_distance_bp=int(l.distance_to_math_candidate_bp),
                     repeat_group_n=len(rr),repeat_envelope_start0=rep_start,repeat_envelope_end0=rep_end,repeat_to_rt_gap_bp=rep_gap,repeat_side_relative_to_rt=side,
                     myrt_model=f.myrt_model,myrt_full_E=d.full_E,myrt_domain_iE=d.i_E,myrt_model_coverage=model_cov,myrt_query_coverage=query_cov,myrt_alignment_aa=f'{d.ali_from}-{d.ali_to}',
                     yxdd_motif_n=len(yxdd),yxdd_motifs=';'.join(f'{pos}:{mot}' for pos,mot in yxdd) if yxdd else 'NONE',yxdd_in_hmm_alignment=';'.join(f'{pos}:{mot}' for pos,mot in motif_in_aln) if motif_in_aln else 'NONE'))
    if len(rr): reps.append(rr)
    motifs.append(dict(target_id=target,record_key=p.record_key,sequence=aa,yxdd=';'.join(f'{pos}:{mot}' for pos,mot in yxdd) if yxdd else 'NONE'))
A=pd.DataFrame(rows)
A.to_csv(out/'THREE_LOCUS_ARCHITECTURE_RT_AUDIT.tsv',sep='\t',index=False)
if reps: pd.concat(reps).to_csv(out/'THREE_LOCUS_REPEAT_GROUPS.tsv',sep='\t',index=False)
pd.DataFrame(motifs).to_csv(out/'THREE_LOCUS_RT_MOTIF_AUDIT.tsv',sep='\t',index=False)
# pair identity and architecture offset for NC_031007 vs NC_048173
p1=A[A.accession=='NC_031007.1'].iloc[0]; p2=A[A.accession=='NC_048173.1'].iloc[0]
pair={
 'same_broad_rt_target': bool(F[F.accession=='NC_031007.1'].broad_rt_target.iloc[0]==F[F.accession=='NC_048173.1'].broad_rt_target.iloc[0]),
 'rt_coordinate_offset_start_bp': int(p2.rt_start0-p1.rt_start0),
 'rt_coordinate_offset_end_bp': int(p2.rt_end0-p1.rt_end0),
 'repeat_envelope_offset_start_bp': int(p2.repeat_envelope_start0-p1.repeat_envelope_start0),
 'repeat_envelope_offset_end_bp': int(p2.repeat_envelope_end0-p1.repeat_envelope_end0),
 'interpretation':'Two genomic instances of one highly conserved architecture lead; do not count as two independent protein-family discoveries.'
}
(out/'PAIR_A_REPLICATION_AUDIT.json').write_text(json.dumps(pair,indent=2)+'\n')
# concise text report
lines=['# Stage05D2 preliminary three-locus architecture + RT-validity audit','',
'Discovery boundary: all three loci were generated by the frozen DNA-only ACNV + generic recurrence detector. Protein/HMM/motif evidence below is post-discovery interpretation and does not modify mathematical scores or ranks.','']
for _,r in A.iterrows():
    lines += [f"## {r.candidate_key}",
              f"- Math rank {r.region_rank}; region [{r.candidate_start0},{r.candidate_end0}); stratum {r.stratum}.",
              f"- Broad-RT-like ORF {r.rt_record_key}: [{r.rt_start0},{r.rt_end0}) {r.rt_strand} strand, {r.rt_aa_len} aa; model {r.myrt_model}; full E={r.myrt_full_E:.3g}.",
              f"- HMM model coverage={r.myrt_model_coverage:.3f}; query coverage={r.myrt_query_coverage:.3f}; aligned aa {r.myrt_alignment_aa}.",
              f"- Prediction-relative noncoding repeat envelope [{r.repeat_envelope_start0},{r.repeat_envelope_end0}); gap to ORF={r.repeat_to_rt_gap_bp} bp; strand-aware side={r.repeat_side_relative_to_rt}.",
              f"- YxDD scan: {r.yxdd_motifs}; motif inside HMM-aligned segment: {r.yxdd_in_hmm_alignment}.",
              '- Interpretation: broad RT-like support only; not confirmed RT activity and not new ART identity.','']
(out/'STAGE05D2_PRELIMINARY_REPORT.md').write_text('\n'.join(lines)+'\n')
print(A.to_string(index=False))
