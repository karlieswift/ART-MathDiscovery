#!/usr/bin/env python3
"""Independent phage sample accession screening; no unannotated genomes become ART controls."""
import argparse,csv,json,pathlib,re

def read_tsv(path):
 with open(path,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def find_value(row,names):
 normalized={re.sub('[^a-z0-9]','',k.lower()):v for k,v in row.items() if k}
 for n in names:
  k=re.sub('[^a-z0-9]','',n.lower())
  if k in normalized:return (normalized[k] or '').strip()
 return ''

def run(metadata,exposure,out,min_len=150000,limit=250):
 exposure_rows=read_tsv(exposure);acc={r['accession_version'] for r in exposure_rows};base={x.split('.')[0] for x in acc}
 rows=read_tsv(metadata);eligible=[];audit=[];seen=set()
 for r in rows:
  raw=find_value(r,['Accession','accession_version','genome accession','GenBank accession'])
  accession=raw.split()[0].strip()
  ntext=find_value(r,['Genome Length (bp)','length_bp','genome_length','Genome Length'])
  try:n=int(ntext.replace(',',''))
  except ValueError:n=0
  qual=find_value(r,['CheckV quality','checkv_quality','Genome quality'])
  description=find_value(r,['Description','genome_description','phage_name'])
  acbase=accession.split('.')[0];exposed=accession in acc or acbase in base
  if not accession or n<min_len or exposed:continue
  if accession in seen:continue
  seen.add(accession)
  version_ok=bool(re.fullmatch(r'[A-Z]{1,5}_?\d+\.\d+',accession))
  d={'accession_or_base':accession,'genome_length_bp':n,'quality_from_source':qual,'description':description,
     'exact_version_in_source':version_ok,'historical_accession_or_base_overlap':False,
     'status':'PRIORITY_FETCH_AUTHOR_ANNOTATION_NEEDED' if version_ok else 'RESOLVE_INSDC_ACCESSION_VERSION_FIRST',
     'is_ART_positive':'UNKNOWN','is_ART_negative':'UNKNOWN','family_independence':'NOT_YET_TESTED',
     'source_data_file':str(metadata),'metadata_selection_rule':f'genome_length_bp>={min_len}; prior exposure base absent'}
  eligible.append(d)
 eligible.sort(key=lambda d:(-d['genome_length_bp'],d['accession_or_base']))
 selected=eligible[:limit]
 out=pathlib.Path(out);out.mkdir(parents=True,exist_ok=True)
 if selected:
  with open(out/'N1_UNLABELED_PHAGE_PRIORITY_QUEUE.tsv','w',encoding='utf-8',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(selected[0]),delimiter='\t');w.writeheader();w.writerows(selected)
 else:(out/'N1_UNLABELED_PHAGE_PRIORITY_QUEUE.tsv').write_text('accession_or_base\tstatus\n',encoding='utf-8')
 result={'n_metadata_rows':len(rows),'n_nonexposed_candidate_genomes_at_least_min_length':len(eligible),'n_priority_queue':len(selected),'n_independent_ART_truth':0,'gate':'CANDIDATE_GENOMES_ONLY_NO_AUTHOR_TRUTH'}
 (out/'N1_METADATA_SCREEN_STATUS.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--metadata',required=True);p.add_argument('--exposure',required=True);p.add_argument('--out',required=True);p.add_argument('--min-length',type=int,default=150000);p.add_argument('--limit',type=int,default=250)
 a=p.parse_args();print(json.dumps(run(a.metadata,a.exposure,a.out,a.min_length,a.limit),indent=2))
