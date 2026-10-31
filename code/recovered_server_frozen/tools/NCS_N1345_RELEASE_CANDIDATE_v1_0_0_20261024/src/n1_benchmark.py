#!/usr/bin/env python3
"""ART N1: rigorously guarded admission and no-truth-leak evaluator.

The algorithm will NEVER certify data independence by itself. Manual approval
MUST be recorded outside the evaluation prediction file.
"""
import argparse, csv, hashlib, json, pathlib, statistics, gzip
REQUIRED=['source_id','accession_version','genome_fasta','genome_sha256','publication_doi','array_start0','array_end0','rt_start0','rt_end0','annotation_evidence_uri','prior_contact_status','family_independence_status','curation_reviewer','curation_approval']

def read(p):
 with open(p,encoding='utf-8',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
 with open(p,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def file_sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def fasta_sequence(p,accession_version=None):
 lines=[];head=[]
 with (gzip.open(p,'rt',encoding='utf-8') if str(p).endswith('.gz') else open(p,encoding='utf-8')) as f:
  for line in f:
   line=line.strip()
   if line.startswith('>'):head.append(line)
   elif line:lines.append(line.upper())
 if len(head)!=1:raise ValueError('not single-record FASTA')
 if accession_version and head[0].split()[0][1:]!=accession_version:
  raise ValueError('FASTA header accession.version mismatch')
 s=''.join(lines)
 if not s or set(s)-set('ACGTNRYKMSWBDHV'):raise ValueError('FASTA invalid IUPAC DNA')
 return s

def admit(input_file,exposure_file,out):
 exposed={r['accession_version'] for r in read(exposure_file)}
 rows=[];accepted=[]
 for r in read(input_file):
  flags=[f'MISSING:{k}' for k in REQUIRED if not r.get(k,'').strip()]
  acc=r.get('accession_version','').strip()
  if acc in exposed: flags.append('HISTORICAL_EXPOSURE_ACCESSION')
  for k,expected in [('prior_contact_status','CERTIFIED_NEVER_EXPOSED_BEFORE_PRE_FREEZE'),('family_independence_status','CERTIFIED_FAMILY_INDEPENDENT'),('curation_approval','APPROVED_BY_HUMAN_CURATOR')]:
   if r.get(k)!=expected:flags.append('NOT_APPROVED:'+k)
  try:
   a,b,c,d=[int(r[z]) for z in ['array_start0','array_end0','rt_start0','rt_end0']]
   if not(0<=a<b and 0<=c<d):flags.append('INVALID_INTERVAL')
  except (ValueError,KeyError,TypeError):flags.append('INVALID_INTERVAL');a=b=c=d=0
  p=pathlib.Path(r.get('genome_fasta',''))
  if not p.is_file():flags.append('FASTA_NOT_FOUND')
  else:
   if file_sha(p)!=r.get('genome_sha256'):flags.append('FASTA_SHA256_MISMATCH')
   try:
    seq=fasta_sequence(p,acc)
    if max(b,d)>len(seq):flags.append('INTERVAL_EXCEEDS_GENOME')
   except Exception as ex:flags.append('INVALID_FASTA:'+str(ex))
  state='HOLD' if flags else 'READY_FOR_BLIND_EVALUATION_HUMAN_CERTIFIED'
  rows.append({'source_id':r.get('source_id',''),'accession_version':acc,'admission_status':state,'issues':';'.join(flags)})
  if not flags:accepted.append(acc)
  # Duplicate accession is not an independent test unit.
 if len(accepted)!=len(set(accepted)):
  for r in rows:
   if r['accession_version'] in accepted:r['admission_status']='HOLD';r['issues']+=';DUPLICATE_ACCESSION'
  accepted=[]
 write(out,rows,['source_id','accession_version','admission_status','issues'])
 return {'n_source_rows':len(rows),'n_admitted_genomes':len(accepted),'primary_benchmark_complete':False,'status':'HOLD_NO_APPROVED_NEW_GROUND_TRUTH' if not accepted else 'READY_TO_SCORE_BLIND_PREDFILES','note':'No pre-freeze provenance can be proved programmatically.'}

def evaluate(admission_tsv,truth_tsv,preds_tsv,out):
 adm={r['accession_version'] for r in read(admission_tsv) if r['admission_status']=='READY_FOR_BLIND_EVALUATION_HUMAN_CERTIFIED'}
 truth=read(truth_tsv);pred=read(preds_tsv)
 if not adm or any(r['accession_version'] not in adm for r in truth):raise ValueError('No admitted primary data, or truth outside admitted group.')
 methods=sorted({r['method'] for r in pred}); by={}
 for p in pred:
  key=(p['accession_version'],p['method']);a=int(p['start0']);b=int(p['end0']);rank=int(p['rank'])
  if a>=b or a<0 or rank<1:raise ValueError('Invalid prediction coordinates/rank')
  by.setdefault(key,[]).append((rank,a,b))
 outrows=[]
 for t in truth:
  acc=t['accession_version'];src=t['source_id'];a=int(t['array_start0']);b=int(t['array_end0'])
  for method in methods:
   rank=0
   for rr,xx,yy in sorted(by.get((acc,method),[])):
    if xx<b and yy>a:rank=rr;break
   outrows.append({'source_id':src,'accession_version':acc,'method':method,'first_hit_rank':rank,'hit_at_1':int(bool(rank and rank<=1)),'hit_at_3':int(bool(rank and rank<=3)),'hit_at_5':int(bool(rank and rank<=5)),'reciprocal_rank':1.0/rank if rank else 0.0,'n_candidates':len(by.get((acc,method),[]))})
 write(out,outrows,list(outrows[0]))
 return {m:{'n_truth_units':sum(z['method']==m for z in outrows),'hits_at_5':sum(z['hit_at_5'] for z in outrows if z['method']==m),'mrr':statistics.mean(z['reciprocal_rank'] for z in outrows if z['method']==m)} for m in methods}

def main():
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='cmd',required=True)
 a=s.add_parser('admit');a.add_argument('--intake',required=True);a.add_argument('--exposure',required=True);a.add_argument('--out',required=True)
 e=s.add_parser('evaluate');e.add_argument('--admission',required=True);e.add_argument('--truth',required=True);e.add_argument('--predictions',required=True);e.add_argument('--out',required=True)
 x=p.parse_args();r=admit(x.intake,x.exposure,x.out) if x.cmd=='admit' else evaluate(x.admission,x.truth,x.predictions,x.out);print(json.dumps(r,indent=2))
if __name__=='__main__':main()
