#!/usr/bin/env python3
"""US-align 67-query reference-normalization sensitivity, no expert certification."""
import argparse,csv,collections,statistics,math,json,random
from pathlib import Path
from openpyxl import load_workbook

def write(path,rows):
 if not rows:raise ValueError('Empty result')
 with open(path,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def num(v):return float(v) if isinstance(v,(int,float)) else None
def median(x):return round(statistics.median(x),6) if x else ''

def analyze(review_file,pairs_file,out):
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 if pairs_file is None:
  w=load_workbook(review_file,read_only=True,data_only=True)
  v=w['Review69'].iter_rows(values_only=True);h=next(v);review=[dict(zip(h,row)) for row in v if row[0]]
  p=w['Pairs2144'].iter_rows(values_only=True);headers=next(p);pairs=[dict(zip(headers,row)) for row in p if row[0]]
 else:
  with open(review_file,newline='') as f:review=list(csv.DictReader(f,delimiter='\t'))
  with open(pairs_file,newline='') as f:pairs=list(csv.DictReader(f,delimiter='\t'))
  for r in review:
   for k in ['pairs completed','Query length (aa)']:
    if r.get(k):r[k]=int(r[k])
  for r in pairs:
   r['Process exit code']=int(r['Process exit code']) if r.get('Process exit code') else -1
   for k in ['TM score, query normalization','TM score, reference normalization','Minimum normalized TM score','Mean normalized TM score','Aligned fraction of query','Aligned fraction of reference','Fraction of paired query residues with pLDDT >= 70']:
    if r.get(k):r[k]=float(r[k])
 byq=collections.defaultdict(list)
 for r in pairs:
  if r['status']!='COMPLETE_FULL_LENGTH_ORTHOGONAL_ALIGNMENT' or r['Process exit code']!=0:raise ValueError('Non-complete alignment found')
  byq[r['Unique query key']].append(r)
 if len(review)!=69 or len(pairs)!=2144 or len(byq)!=67:raise ValueError('Unexpected 69 or 2144 input counts')
 rows=[]
 for z in review:
  key=z['Unique query key'];group=byq.get(key,[]); done=len(group)==32
  if (z['pairs completed']==32)!=done:raise ValueError('Mismatched pair count '+str(key))
  def getbest(f):return max(group,key=lambda r:(float(r[f]),r['Reference ID'])) if group else None
  best_min=getbest('Minimum normalized TM score')
  best_q=getbest('TM score, query normalization')
  best_ref=getbest('TM score, reference normalization')
  best_mean=getbest('Mean normalized TM score')
  sorted_min=sorted(group,key=lambda r:float(r['Minimum normalized TM score']),reverse=True)
  margin=float(sorted_min[0]['Minimum normalized TM score'])-float(sorted_min[1]['Minimum normalized TM score']) if len(sorted_min)>=2 else None
  flags=[]
  if not done:flags.append('HISTORICAL_OOM_NO_STRUCTURE')
  else:
   if float(best_min['Minimum normalized TM score'])<0.3:flags.append('LOW_GLOBAL_TM_MIN_DESCRIPTIVE')
   if float(best_min['Aligned fraction of query'])<0.5:flags.append('QUERY_COVERAGE_LT50_DESCRIPTIVE')
   if float(best_min['Fraction of paired query residues with pLDDT >= 70'])<0.5:flags.append('PAIRED_PLDDT_CONFIDENCE_LT50_DESCRIPTIVE')
   if len(set([best_min['Reference ID'],best_q['Reference ID'],best_ref['Reference ID']]))>1:flags.append('NORMALIZATION_CHANGES_TOP_REFERENCE')
   if margin is not None and margin<0.02:flags.append('TOP_TWO_TM_MIN_WITHIN_0P02_DESCRIPTIVE')
  rows.append({'membership_id':z['Membership ID'],'cohort':z['Cohort'],'query_id':z['Query ID'],'query_length_aa':z['Query length (aa)'],'archived_tier_unchanged':z['Archived tier'],'pairs_count':len(group),'best_ref_by_tm_min':best_min['Reference ID'] if done else '', 'best_ref_by_tm_query':best_q['Reference ID'] if done else '', 'best_ref_by_tm_reference':best_ref['Reference ID'] if done else '', 'best_ref_by_tm_mean':best_mean['Reference ID'] if done else '', 'tm_min_of_best_pair':num(best_min['Minimum normalized TM score']) if done else '', 'best_query_fraction':num(best_min['Aligned fraction of query']) if done else '', 'best_reference_fraction':num(best_min['Aligned fraction of reference']) if done else '', 'paired_query_plddt_ge70_fraction':num(best_min['Fraction of paired query residues with pLDDT >= 70']) if done else '', 'top_two_tm_min_margin':round(margin,6) if done else '', 'different_normalization_winner':int(len(set([best_min['Reference ID'],best_q['Reference ID'],best_ref['Reference ID']]))>1) if done else '', 'descriptive_flags':';'.join(flags) if flags else 'NO_FLAG','expert_review':'PENDING_HUMAN_REVIEW','biochemical_function':'NOT_TESTED','scope':'NO_NEW_FOLDSEEK_TIER_CLASSIFICATION'})
 write(out/'N4_QUERY_NORMALIZATION_AUDIT_69.tsv',rows)
 group_summary=[]
 for cohort in sorted({r['cohort'] for r in rows}):
  x=[r for r in rows if r['cohort']==cohort];k=[r for r in x if r['pairs_count']==32]
  group_summary.append({'cohort':cohort,'n_memberships':len(x),'n_full_pair_queries':len(k),'n_oom':len(x)-len(k),'median_tm_min':median([r['tm_min_of_best_pair'] for r in k]),'median_best_aligned_query_fraction':median([r['best_query_fraction'] for r in k]),'median_two_best_tm_min_gap':median([r['top_two_tm_min_margin'] for r in k]),'n_normalization_changes_best_reference':sum(r['different_normalization_winner'] for r in k),'n_low_tm_min_lt0p3':sum(r['tm_min_of_best_pair']<0.3 for r in k),'n_low_query_coverage_lt0p5':sum(r['best_query_fraction']<0.5 for r in k),'n_low_plddt_paired_lt0p5':sum(r['paired_query_plddt_ge70_fraction']<0.5 for r in k),'scope':'Descriptive, not independent replication or biological activity'})
 write(out/'N4_GROUP_SUMMARY.tsv',group_summary)
 # Four archived strict positives previously prioritised by full-length alignment flags.
 positives=[r for r in rows if r['cohort']=='strict_mirror_54Gb' and r['archived_tier_unchanged'] in {'STRONG','SUGGESTIVE'}]
 flagged=[r for r in positives if r['tm_min_of_best_pair']!='' and (r['tm_min_of_best_pair']<0.3 or r['best_query_fraction']<0.5 or r['paired_query_plddt_ge70_fraction']<0.5)]
 write(out/'N4_STRICT_24_PRIORITIZED.tsv',[dict(r) for r in positives])
 report={'gate':'PASS_N4_67x32_INPUT_CONCORDANCE_AND_NORMALIZATION_AUDIT','n_memberships':69,'n_pairs':2144,'n_structures_available':67,'n_oom':2,'strict_positive_count':len(positives),'strict_prioritized_any_original_descriptive_flag':len(flagged),'has_human_expert_curation':False,'tier_reclassification':False,'normalization_winner_changes':sum(r['different_normalization_winner'] for r in rows if isinstance(r['different_normalization_winner'],int))}
 (out/'N4_STATUS.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--review',required=True);p.add_argument('--pairs');p.add_argument('--out',required=True);a=p.parse_args();print(json.dumps(analyze(a.review,a.pairs,a.out),indent=2))
