#!/usr/bin/env python3
"""Post-score focused comparison and score-blind RNA input preparation."""
from pathlib import Path
import argparse,itertools,json
from common import *
from Bio import Align,__version__ as bio_version
from Bio.Align import substitution_matrices

ORDER=['OR836606.1','LC680885.1','MZ779063.1','MN091626.1','OR455461.1','NC_016163.1','MT006214.1','OR770614.1','NC_067211.1','OK490448.1']
CLOSE=set(ORDER[:5]+['OR770614.1'])

def align_metrics(x,y,protein):
 a=Align.PairwiseAligner();a.mode='global'
 if protein:a.substitution_matrix=substitution_matrices.load('BLOSUM62');a.open_gap_score=-11;a.extend_gap_score=-1
 else:a.match_score=2;a.mismatch_score=-3;a.open_gap_score=-5;a.extend_gap_score=-2
 z=a.align(x,y)[0];c=z.coordinates;matches=mm=gap=gapruns=paired=0;variants=[]
 for i in range(c.shape[1]-1):
  p,q=map(int,c[:,i]);u,v=map(int,c[:,i+1]);dx=u-p;dy=v-q
  if dx and dy:
   if dx!=dy:raise ValueError('unexpected diagonal')
   paired+=dx
   for j,(xx,yy) in enumerate(zip(x[p:u],y[q:v])):
    if xx==yy:matches+=1
    else:
     mm+=1;variants.append({'type':'SUBSTITUTION','a_start':p+j,'a_end':p+j+1,'b_start':q+j,'b_end':q+j+1,'a_seq':xx,'b_seq':yy})
  else:
   gap+=max(dx,dy);gapruns+=1
   variants.append({'type':'GAP_IN_B' if dx else 'GAP_IN_A','a_start':p,'a_end':u,'b_start':q,'b_end':v,'a_seq':x[p:u],'b_seq':y[q:v]})
 cols=paired+gap
 return {'a_length':len(x),'b_length':len(y),'alignment_score':float(z.score),'matching_columns':matches,'mismatch_columns':mm,'gap_columns':gap,'gap_runs':gapruns,'identity_all_columns_pct':100*matches/cols,'paired_a_fraction':paired/len(x),'paired_b_fraction':paired/len(y),'exact_sequence_identical':x==y},variants

def prepare(root,out):
 root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
 E=root/'inputs/evidence';F=root/'inputs/focused_sequences'
 policy=read_json(root/'config/POLICY.json')
 co={r['sequence_id']:r for r in read_tsv(E/'TEN_LOCUS_INTEGRATED_COMPARISON.tsv')}
 loci={r['sequence_id']:r for r in read_tsv(E/'CORE10_PROFILE_DETAIL.tsv')}
 if set(co)!=set(ORDER) or set(loci)!=set(ORDER):raise ValueError('ten-locus set changed')
 up=fasta(E/'TEN_RT_ORIENTED_UPSTREAM_INTERGENIC.fna');upmeta={r['sequence_id']:r for r in read_tsv(E/'TEN_RT_ORIENTED_UPSTREAM_INTERVALS.tsv')}
 proteins=fasta(F/'TEN_RT_PROTEINS.faa');neighbors=fasta(F/'ALL_TEN_LOCI_NEIGHBOR_ORFS.faa')
 gs=read_tsv(E/'FOCUSED_NEIGHBOR_ORF_REGISTRY.tsv')
 for g in gs:
  s=neighbors[g['record_key']]
  if len(s)!=int(g['aa_length']) or sha_text(s)!=g['aa_sha256']:raise ValueError('neighbor source mismatch')
 contexts={};bounds={};seen_genes={}
 for m in read_tsv(E/'TEN_RT_CONTEXT_SEQUENCE_MANIFEST.tsv'):
  p=F/m['file'];d=fasta(p)
  if len(d)!=1:raise ValueError('context must be single record')
  s=next(iter(d.values()));acc=m['sequence_id'];a,b=int(m['start0']),int(m['end0'])
  if sha_file(p)!=m['fasta_sha256'] or sha_text(s)!=m['sequence_sha256'] or len(s)!=b-a or set(s)-set('ACGT'):raise ValueError('context mismatch '+acc)
  contexts[acc]=s;bounds[acc]=(a,b)
  l=loci[acc];u=upmeta[acc];nt=s[int(u['start0'])-a:int(u['end0'])-a]
  if u['rt_strand']=='-':nt=rc(nt)
  if nt!=up[acc] or sha_text(nt)!=u['sha256']:raise ValueError('upstream exact replay failed '+acc)
  seen_genes[acc]=[(int(g['start0']),int(g['end0'])) for g in gs if g['sequence_id']==acc]+[(int(l['start0']),int(l['end0']))]
 ledger=[]
 for acc in ORDER:
  r=co[acc];l=loci[acc];inside=bool(r['math_candidate_key'])
  group=('A_CAPTURED_CLOSE_ART_REFERENCE_RELATIVE' if acc in CLOSE and inside else 'D_NOT_CAPTURED_CLOSE_ART_REFERENCE_RELATIVE' if acc in CLOSE else 'B_CAPTURED_DISTANT_COOCCURRENCE' if acc=='NC_016163.1' else 'C_CAPTURED_RT_PROFILE_WITHOUT_QUALIFYING_ARRAY' if inside else 'E_NOT_CAPTURED_RT_PROFILE_OTHER')
  ledger.append({'sequence_id':acc,'whole_genome_RT_core_GA_locus':1,'math_neighborhood_captured':inside,'math_candidate_key':r['math_candidate_key'],'math_region_rank':r['math_region_rank'],'C1_exact_plus_RT_label':acc in set(ORDER[:4]+['NC_016163.1']),'C1_approximate_plus_RT_label':acc=='OR455461.1','C2A_close_three_component_reference_relative':acc in CLOSE,'review_group':group,'rt_start0':l['start0'],'rt_end0':l['end0'],'rt_strand':l['strand'],'rt_aa_length':l['aa_length'],'family_novelty_confirmed':False,'wet_lab_ready':False})
 write_tsv(out/'TEN_LOCUS_UNAMBIGUOUS_LEDGER.tsv',ledger)
 counters={'whole_genome_RT_core_GA_loci':10,'math_neighborhood_captured':7,'math_neighborhood_not_captured':3,'C1_exact_joint_labels':5,'C1_approximate_joint_labels':1,'C2A_captured_close_reference_relatives':5,'C2A_captured_distant_cooccurrence':1,'C2A_captured_RT_only':1,'C2A_uncaptured_close_reference_relatives':1,'all_close_reference_relatives':6,'confirmed_new_ART_families':0}
 write_json(out/'COUNT_DEFINITIONS.json',counters)
 # New continuous comparisons: all pairs of the SIX already reviewed close loci.
 # This selection is post-hoc, descriptive, not an accuracy or discovery test.
 comps=[];changes=[];seqout={}
 for i,j in itertools.combinations([x for x in ORDER if x in CLOSE],2):
  for kind in ['RT_PROTEIN','PARTNER_LIKE_PROTEIN','FULL_UPSTREAM_DNA']:
   if kind=='RT_PROTEIN':s,t=proteins[co[i]['rt_record_key']],proteins[co[j]['rt_record_key']]
   elif kind=='PARTNER_LIKE_PROTEIN':s,t=neighbors[co[i]['best_neighbor_ORF']],neighbors[co[j]['best_neighbor_ORF']]
   else:s,t=up[i],up[j]
   ms,vs=align_metrics(s,t,kind!='FULL_UPSTREAM_DNA')
   comps.append(dict(a=i,b=j,component=kind,**ms))
   if kind=='FULL_UPSTREAM_DNA':
    for v in vs:changes.append(dict(a=i,b=j,coordinate_system='zero_based_RT_oriented_upstream_not_genome',**v))
   seqout[i+'__'+kind]=s;seqout[j+'__'+kind]=t
 write_tsv(out/'SIX_CLOSE_LOCUS_GLOBAL_COMPONENT_COMPARISON.tsv',comps)
 write_tsv(out/'SIX_CLOSE_UPSTREAM_VARIANT_COLUMNS.tsv',changes)
 uniq=[]
 for k in ['RT_PROTEIN','PARTNER_LIKE_PROTEIN','FULL_UPSTREAM_DNA']:
  sd={acc:seqout[acc+'__'+k] for acc in CLOSE}
  uniq.append({'component':k,'loci':6,'unique_exact_sequences':len(set(sd.values())),'independent_families_inferred':False})
 write_tsv(out/'EXACT_SEQUENCE_DIVERSITY.tsv',uniq)
 # All 10 fixed RT-centered intergenic contexts. No folding score is read here.
 occurrences={};targets={};excluded=[];ic=[]
 def add(acc,a,b,kind,nt=None,reference=False,tag=''):
  if b-a<60:
   excluded.append({'sequence_id':acc,'start0':a,'end0':b,'kind':kind,'reason':'SHORTER_THAN_60NT_NOT_A_BIOLOGICAL_NEGATIVE'});return
  if nt is None:
   ca,cb=bounds[acc]
   if not ca<=a<b<=cb:raise ValueError('out of fixed context '+acc)
   nt=contexts[acc][a-ca:b-ca]
  if len(nt)!=b-a or set(nt)-set('ACGT'):raise ValueError('bad RNA interval source')
  for strand,s in [('+',nt),('-',rc(nt))]:
   rna=s.replace('T','U');tid='RNA_'+sha_text(rna);targets[tid]=rna;key=(acc,a,b,strand)
   if key not in occurrences:
    occurrences[key]={'sequence_id':acc,'start0':a,'end0':b,'genomic_strand_hypothesis':strand,'relative_to_RT':'REGISTERED_REFERENCE' if reference else 'SAME_AS_RT' if strand==loci[acc]['strand'] else 'OPPOSITE_TO_RT','kinds':set(),'repeat_description_ids':set(),'target_id':tid,'length_nt':b-a,'RNA_sha256':sha_text(rna),'math_neighborhood_captured':False if reference else bool(co[acc]['math_candidate_key']),'registered_reference':reference,'transcript_boundary_proven':False}
   occurrences[key]['kinds'].add(kind)
   if tag:occurrences[key]['repeat_description_ids'].add(tag)
 for acc in ORDER:
  a,b=bounds[acc]
  gaps=complement(a,b,seen_genes[acc])
  for x,y in gaps:
   ic.append({'sequence_id':acc,'start0':x,'end0':y,'length_nt':y-x,'scope':'relative_to_frozen_Prodigal_CDS_not_proven_noncoding_RNA'})
   for p,q,k in windows(x,y):add(acc,p,q,k)
  u=upmeta[acc];a0,b0=int(u['start0']),int(u['end0'])
  if b0-a0<=2048:add(acc,a0,b0,'FULL_IMMEDIATE_UPSTREAM')
 reps=read_tsv(E/'CORE_CANDIDATE_QUALIFYING_REPEAT_DESCRIPTIONS.tsv')
 for r in reps:
  acc=r['sequence_id'];a,b=bounds[acc];iv=json.loads(r['copy_intervals_json'])
  for (p,q),(x,y) in zip(iv,iv[1:]):
   for s,e,k in [(p,x,'REPEAT_PLUS_SPACER'),(q,y,'SPACER_PLUS_REPEAT')]:
    if a<=s<e<=b and e-s<=2048:add(acc,s,e,k,tag=r['repeat_group_id'])
    else:excluded.append({'sequence_id':acc,'start0':s,'end0':e,'kind':k,'reason':'OUTSIDE_FIXED_RT_CONTEXT_OR_TOO_LONG_NOT_A_NEGATIVE'})
 refcore=fasta(E/'REGISTERED_ARRAY_CORE_REFERENCES.fna')
 refm=read_tsv(root/'inputs/reference_sources/REPEAT_REFERENCE_CONTEXT.tsv')
 for rr in refm:
  acc=rr['reference_id'];nt=refcore[acc]
  # Field names are validated against actual registered source below.
  a=int(rr['registered_core_start0']);b=int(rr['registered_core_end0'])
  add(acc,a,b,'FULL_REGISTERED_CORE',nt,True)
  for s,e,k in windows(a,b):add(acc,s,e,'REFERENCE_'+k,nt[s-a:e-a],True)
 rows=[]
 for key in sorted(occurrences):
  r=occurrences[key];r['kinds']=';'.join(sorted(r['kinds']));r['repeat_description_ids']=';'.join(sorted(r['repeat_description_ids']));r={'occurrence_id':'O%06d'%(len(rows)+1),**r};rows.append(r)
 write_tsv(out/'RNA_OCCURRENCE_REGISTRY.tsv',rows)
 write_fasta(out/'RNA_TARGETS.fasta',dict(sorted(targets.items())))
 write_tsv(out/'INTERGENIC_INTERVALS.tsv',ic)
 write_tsv(out/'RNA_SCOPE_EXCLUSIONS.tsv',excluded,['sequence_id','start0','end0','kind','reason'])
 freq=[{'sequence_id':acc,'RNA_occurrences':sum(r['sequence_id']==acc for r in rows),'RNA_unique_targets':len({r['target_id'] for r in rows if r['sequence_id']==acc}),'math_captured':bool(co.get(acc,{}).get('math_candidate_key','')),'reference':acc.startswith('REF_')} for acc in ORDER+list(refcore)]
 write_tsv(out/'RNA_INPUT_COUNTS.tsv',freq)
 freeze={'stage':'Stage05C2B','gate':'PASS_TEN_LOCUS_LEDGER_VARIATION_AND_RNA_INPUT_FREEZE','n_loci':10,'n_references':3,'n_pairwise_component_comparisons':len(comps),'n_RNA_occurrences':len(rows),'n_unique_RNA_targets':len(targets),'folds_with_99_nulls':len(targets)*100,'RNAfold_executed':False,'full_author_catalogue_compared':False,'RNA_targets_sha256':sha_file(out/'RNA_TARGETS.fasta'),'RNA_occurrences_sha256':sha_file(out/'RNA_OCCURRENCE_REGISTRY.tsv'),'policy_sha256':sha_file(root/'config/POLICY.json'),'Bio_version':bio_version,'both_directions_retained':True,'UG27_touched':False,'novel_ART_confirmed':0}
 write_json(out/'INPUT_FREEZE.json',freeze)
 print(json.dumps(freeze,indent=2))
 return freeze
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--root',required=True);a.add_argument('--out',required=True);x=a.parse_args();prepare(x.root,x.out)
