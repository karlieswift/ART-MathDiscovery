#!/usr/bin/env python3
"""Post-freeze 6-genome region/rank stability analysis (not accuracy)."""
import argparse,csv,gzip,collections,random,json,statistics
from pathlib import Path

def read(p):
 f=gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)
 with f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows):
 if not rows:raise ValueError('No rows')
 with open(path,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def union(records):
 merged=[]
 for a,b in sorted(records):
  if a<0 or b<=a:raise ValueError('invalid interval')
  if not merged or a>merged[-1][1]:merged.append([a,b])
  else:merged[-1][1]=max(b,merged[-1][1])
 return merged
def len_union(vals):return sum(b-a for a,b in vals)
def intersects(x,y):return x[0]<y[1] and y[0]<x[1]
def overlap(a,b):
 i=j=s=0
 while i<len(a) and j<len(b):
  s+=max(0,min(a[i][1],b[j][1])-max(a[i][0],b[j][0]))
  if a[i][1]<=b[j][1]:i+=1
  else:j+=1
 return s

def calculate(input_path,out_dir):
 out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
 grouped=collections.defaultdict(list)
 for row in read(input_path):
  grouped[(row['case_id'],float(row['alpha']))].append({'rank':int(row['rank']),'start':int(row['start0']),'end':int(row['end0'])})
 cases=sorted(set(k[0] for k in grouped));alphas=(0.0025,0.005,0.01)
 if len(cases)!=6 or len(grouped)!=18:raise ValueError('Expected exactly six cases x three thresholds')
 results=[];matches=[]
 for case in cases:
  frozen=grouped[(case,0.005)];base=union([(x['start'],x['end']) for x in frozen]);base_len=len_union(base)
  for alpha in alphas:
   test=grouped[(case,alpha)]
   if sorted(x['rank'] for x in test)!=list(range(1,len(test)+1)):raise ValueError(f'invalid ranks {case} {alpha}')
   cur=union([(x['start'],x['end']) for x in test]);tlen=len_union(cur);ol=overlap(base,cur)
   base_top=sorted(frozen,key=lambda z:z['rank'])[:5]
   test_top=sorted(test,key=lambda z:z['rank'])[:5]
   top_overlap=sum(any(intersects((b['start'],b['end']),(t['start'],t['end'])) for t in test_top) for b in base_top)
   exact_region_matching=sum((x['start'],x['end']) in {(y['start'],y['end']) for y in test} for x in frozen)
   retained_top10=[]
   for b in sorted(frozen,key=lambda z:z['rank'])[:10]:
    relevant=[z for z in test if intersects((b['start'],b['end']),(z['start'],z['end']))]
    if relevant:retained_top10.append(abs(b['rank']-min(z['rank'] for z in relevant)))
    matches.append({'case_id':case,'q_threshold':alpha,'frozen_rank':b['rank'],'frozen_start0':b['start'],'frozen_end0':b['end'],'overlap_any_test':int(bool(relevant)),'best_overlapping_test_rank':min((z['rank'] for z in relevant),default=''),'absolute_rank_shift_if_matched':min((abs(b['rank']-z['rank']) for z in relevant),default='')})
   results.append({'case_id':case,'q_threshold':alpha,'n_regions':len(test),'frozen_regions':len(frozen),'frozen_covered_bp':base_len,'test_covered_bp':tlen,'bp_overlap':ol,'bp_recall_from_frozen':ol/base_len,'interval_jaccard':ol/(base_len+tlen-ol),'retained_exact_region_count':exact_region_matching,'frozen_top5_hit_test_top5':top_overlap,'frozen_top10_overlap_any_test':len(retained_top10),'median_abs_rank_shift_top10_matched':statistics.median(retained_top10) if retained_top10 else '', 'scope':'POSTFREEZE_DIAGNOSTIC_NOT_PERFORMANCE'})
 write(out_dir/'N3_ROBUSTNESS_CASES.tsv',results);write(out_dir/'N3_TOP10_MATCHING.tsv',matches)
 summary=[];rgen=random.Random(20261024)
 for a in alphas:
  rr=[r for r in results if r['q_threshold']==a]; vals=[r['interval_jaccard'] for r in rr]
  samples=[statistics.mean([vals[rgen.randrange(len(vals))] for _ in vals]) for j in range(5000)]
  samples.sort()
  summary.append({'q_threshold':a,'n_genomes':len(rr),'total_regions':sum(r['n_regions'] for r in rr),'median_n_regions':statistics.median(r['n_regions'] for r in rr),'mean_interval_jaccard':statistics.mean(vals),'bootstrap_95_lower_descriptive':samples[125],'bootstrap_95_upper_descriptive':samples[4875],'median_top5_retained':statistics.median(r['frozen_top5_hit_test_top5'] for r in rr),'scope':'six previously used genomes; bootstrap across genomes, no new blind truth'})
 write(out_dir/'N3_ROBUSTNESS_SUMMARY.tsv',summary)
 out={'gate':'PASS_N3_POSTFREEZE_RECOMBINATION_AND_RANK_STABILITY','n_genomes':6,'n_combinations':18,'all_identical_at_frozen_threshold':all(r['interval_jaccard']==1 and r['n_regions']==r['frozen_regions'] for r in results if r['q_threshold']==0.005),'bootstrapping':'descriptive n=6, seed 20261024','independent_benchmark':False}
 (out_dir/'N3_STATUS.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--region-ranks',required=True);p.add_argument('--out',required=True);a=p.parse_args();print(json.dumps(calculate(a.region_ranks,a.out),indent=2))
