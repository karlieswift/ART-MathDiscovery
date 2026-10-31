#!/usr/bin/env python3
import argparse,csv,hashlib,json,sys
from pathlib import Path

def read_tsv(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,cols):
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def fasta(p):
 h=None;seq=[]
 with open(p) as f:
  for line in f:
   line=line.strip()
   if not line:continue
   if line.startswith('>'):
    if h is not None:yield h,''.join(seq).rstrip('*').upper()
    h=line[1:].split()[0];seq=[]
   else:seq.append(line)
  if h is not None:yield h,''.join(seq).rstrip('*').upper()
def shaseq(s):return hashlib.sha256(s.encode()).hexdigest()
def one(P,t):
 ds=sorted(d for d in P.iterdir() if d.is_dir() and t in d.name)
 if not ds:raise RuntimeError('cannot locate '+t)
 return ds[-1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--project-root',required=True);ap.add_argument('--home',required=True);a=ap.parse_args();P=Path(a.project_root).resolve();H=Path(a.home).resolve();O=H/'outputs';I=H/'inputs';O.mkdir(exist_ok=True);I.mkdir(exist_ok=True)
 B8=one(P,'STAGE04B8R1_CROSS_ROUTE_RECONCILIATION_AND_UNRESOLVED5_AUDIT');D=one(P,'STAGE05F3V2_60G_DELTA_DISTRIBUTED_ACQUISITION_R2_v1_0_2_20261013')
 tab=B8/'outputs/UNRESOLVED5_ORIGINAL_CONTIG_STOP_BOUNDARY_AUDIT.tsv';rows=read_tsv(tab)
 if len(rows)!=5:raise RuntimeError(f'expected 5 B8R1 rows got {len(rows)}')
 for r in rows:
  if r['original_contig_status']!='ORIGINAL_CONTIG_STOP_BOUNDED_BOTH' or str(r['left_stop_bounded'])!='True' or str(r['right_stop_bounded'])!='True' or int(r['left_extension_aa'])!=0 or int(r['right_extension_aa'])!=0 or int(r['original_contig_recovered_aa_length'])!=int(r['existing_aa_length']):raise RuntimeError('B8R1 completeness condition not met '+r['candidate_id'])
 faa=D/'outputs/F3V2_FORMAL_ALL.faa';need={r['candidate_id'] for r in rows};seq={q:s for q,s in fasta(faa) if q in need}
 if set(seq)!=need:raise RuntimeError(f'formal FASTA missing B8R1 candidates {need-set(seq)}')
 out=[];yd=I/'arch5_boltz_yaml';yd.mkdir(parents=True,exist_ok=True)
 for i,r in enumerate(sorted(rows,key=lambda x:x['delta_operational_family_id']),1):
  q=f'A5Q{i:04d}';s=seq[r['candidate_id']]
  if len(s)!=int(r['existing_aa_length']):raise RuntimeError('AA length mismatch '+r['candidate_id'])
  y=yd/f'{q}.yaml';y.write_text('version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: '+s+'\n')
  out.append({'query_id':q,'branch':'ARCH5_COMPLETE_ORF_EXTENSION','delta_operational_family_id':r['delta_operational_family_id'],'candidate_id':r['candidate_id'],'accession':r['accession'],'aa_length':len(s),'sequence_sha256':shaseq(s),'completeness_evidence':'ORIGINAL_CONTIG_STOP_BOUNDED_BOTH_ZERO_EXTENSION','b9_role':'POST_B8_TECHNICAL_COMPLETE_ORF_STRUCTURE_EXTENSION','yaml_path':str(y.resolve())})
 write_tsv(O/'ARCH5_STRUCTURE_MANIFEST.tsv',out,list(out[0]))
 gate={'stage':'Stage05F3V2-Stage04B9B','gate':'PASS_STAGE04B9B_ARCH5_COMPLETE_ORF_STRUCTURE_INPUT_PREPARED','gate_pass':True,'target_n':5,'all_original_contig_stop_bounded_both':True,'all_zero_extension':True,'selection_changed':False,'retuning_performed':False,'structure_run_performed':False}
 with open(O/'STAGE04B9B_PREPARE_GATE.json','w') as f:json.dump(gate,f,indent=2);f.write('\n')
 print(json.dumps(gate,indent=2));print('PASS_STAGE04B9B_ARCH5_COMPLETE_ORF_STRUCTURE_INPUT_PREPARED')
if __name__=='__main__':
 try:main()
 except Exception as e:print('[HOLD]',e,file=sys.stderr);raise
