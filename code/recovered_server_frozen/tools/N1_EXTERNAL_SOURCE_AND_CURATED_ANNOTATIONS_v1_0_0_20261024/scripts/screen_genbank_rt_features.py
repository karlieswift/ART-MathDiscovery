#!/usr/bin/env python3
"""Screen author GenBank RT CDS annotations; outputs are candidates only, never ART truth.

Biopython is optional for this additional step: pip install biopython.
"""
import argparse,csv,gzip,hashlib,json,pathlib,re

def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def run(gbdir,exposure,out,max_records=500):
 try:
  from Bio import SeqIO
 except ImportError:
  raise SystemExit('Biopython needed for GBFF screening: python -m pip install biopython')
 with open(exposure,encoding='utf-8') as f:
  data=list(csv.DictReader(f,delimiter='\t'))
 exposure_acc={r['accession_version'] for r in data};exp_base={x.split('.')[0] for x in exposure_acc}
 files=[p for p in pathlib.Path(gbdir).rglob('*') if p.is_file() and p.name.lower().endswith(('.gbff','.gbk','.gb','.gbff.gz','.gbk.gz','.gb.gz'))]
 files=sorted(files);out=pathlib.Path(out);out.mkdir(parents=True,exist_ok=True)
 rows=[];count=0
 for file in files:
  op=gzip.open if file.name.endswith('.gz') else open
  with op(file,'rt') as f:
   for rec in SeqIO.parse(f,'genbank'):
    count+=1
    if count>max_records:raise RuntimeError('Maximum GenBank record count exceeded; narrow input directory')
    acc=rec.id;acc_base=acc.split('.')[0]
    exposed=acc in exposure_acc or acc_base in exp_base
    for ft in rec.features:
     if ft.type!='CDS':continue
     prod=' | '.join(ft.qualifiers.get('product',[])+ft.qualifiers.get('note',[]))
     if not re.search(r'(reverse transcriptase|RNA.?directed DNA polymerase|retron)',prod,re.I):continue
     start=int(ft.location.start);end=int(ft.location.end);strand=ft.location.strand
     rows.append({'accession_version':acc,'genbank_source_file':str(file),'source_sha256':sha(file),
                  'genome_length_bp':len(rec.seq),'feature_start0':start,'feature_end0':end,'strand':strand,
                  'annotated_product':prod[:350],'protein_id':';'.join(ft.qualifiers.get('protein_id',[])),
                  'historically_exposed':exposed,'primary_ART_status':'UNANNOTATED_RT_ONLY_NOT_ART_TRUTH',
                  'required_next_step':'Independent manual noncoding array+partner coordinates and family exposure audit'})
 fields=['accession_version','genbank_source_file','source_sha256','genome_length_bp','feature_start0','feature_end0','strand','annotated_product','protein_id','historically_exposed','primary_ART_status','required_next_step']
 with open(out/'N1_NEW_GENBANK_RT_CDS_CANDIDATES.tsv','w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
 ret={'genbank_files_found':len(files),'genbank_records_reviewed':count,'rt_cds_feature_candidates':len(rows),'new_independent_art_truth':0,'status':'RT_CDS_SCREENED_PENDING_ARRAY_ANNOTATION'}
 (out/'N1_GENBANK_SCREEN_STATUS.json').write_text(json.dumps(ret,indent=2)+'\n',encoding='utf-8')
 return ret
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--gbff-dir',required=True);p.add_argument('--exposure',required=True);p.add_argument('--out',required=True);p.add_argument('--max-records',type=int,default=500)
 a=p.parse_args();print(json.dumps(run(a.gbff_dir,a.exposure,a.out,a.max_records),indent=2))
