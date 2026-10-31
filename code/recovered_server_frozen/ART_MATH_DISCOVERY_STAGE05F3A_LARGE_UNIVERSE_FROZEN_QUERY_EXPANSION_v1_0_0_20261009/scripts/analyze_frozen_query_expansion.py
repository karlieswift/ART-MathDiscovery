#!/usr/bin/env python3
import argparse,csv,hashlib,json,math,re,subprocess,sys
from collections import defaultdict
from pathlib import Path

FORMAL=1e-4;EXPLORE=1e-2
RC=str.maketrans('ACGTN','TGCAN')
CODON={
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

def sha256_file(p):
 h=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def read_fasta(p):
 d={};h=None;s=[]
 with open(p) as f:
  for line in f:
   line=line.strip()
   if not line:continue
   if line.startswith('>'):
    if h is not None:d[h]=''.join(s)
    h=line[1:].split()[0];s=[]
   else:s.append(line.upper())
 if h is not None:d[h]=''.join(s)
 return d

def tr(nt):return ''.join(CODON.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))

def six_segments(seq,minaa=45):
 L=len(seq)
 for strand in ('+','-'):
  nt=seq if strand=='+' else seq.translate(RC)[::-1]
  for frame in (0,1,2):
   aa=tr(nt[frame:]);aa0=0
   for part in aa.split('*'):
    aa1=aa0+len(part)
    if len(part)>=minaa:
     n0=frame+aa0*3;n1=frame+aa1*3
     if strand=='+':s,e=n0,n1
     else:s,e=L-n1,L-n0
     yield strand,frame,s,e,part
    aa0=aa1+1

def parse_m8(p):
 fields=['query','target','evalue','bits','qcov','alnlen','qlen','tlen','fident','qstart','qend','tstart','tend','taln']
 with open(p) as f:
  for r in csv.DictReader(f,delimiter='\t',fieldnames=fields):yield r

def parse_domtbl(p):
 out=[]
 if not Path(p).exists():return out
 with open(p) as f:
  for line in f:
   if not line.strip() or line.startswith('#'):continue
   x=line.rstrip().split(maxsplit=22)
   if len(x)<23:continue
   out.append({'hmm':x[0],'hmm_len':int(x[2]),'query':x[3],'query_len':int(x[5]),'full_e':float(x[6]),'i_e':float(x[12]),'score':float(x[13]),'hmm_from':int(x[15]),'hmm_to':int(x[16]),'ali_from':int(x[17]),'ali_to':int(x[18]),'acc':float(x[21])})
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--queries',required=True);ap.add_argument('--query-table',required=True);ap.add_argument('--components',required=True);ap.add_argument('--external-fasta',required=True);ap.add_argument('--m8',required=True);ap.add_argument('--hmm',required=True);ap.add_argument('--threads',type=int,default=32);ap.add_argument('--outdir',required=True)
 a=ap.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
 qseq=read_fasta(a.queries);qt={r['query_id']:r for r in csv.DictReader(open(a.query_table),delimiter='\t')};
 if len(qseq)!=30 or set(qseq)!=set(qt):raise RuntimeError('query lock mismatch')
 comp={};
 for r in csv.DictReader(open(a.components),delimiter='\t'):
  for q in r['query_ids'].split(';'):comp[q]=r['f2d_component_id']
 if len(set(comp.values()))!=26 or set(comp)!=set(qseq):raise RuntimeError('component lock mismatch')
 rows=list(parse_m8(a.m8));hit_targets={r['target'] for r in rows};contigs={};h=None;s=[]
 with open(a.external_fasta) as f:
  for line in f:
   line=line.strip()
   if not line:continue
   if line.startswith('>'):
    if h is not None and h in hit_targets:contigs[h]=''.join(s)
    h=line[1:].split()[0];s=[]
   else:
    if h in hit_targets:s.append(line.upper())
  if h is not None and h in hit_targets:contigs[h]=''.join(s)
 missing=hit_targets-set(contigs)
 if missing:raise RuntimeError(f'missing {len(missing)} hit target contigs; first={sorted(missing)[:5]}')
 recovered={};annot=[]
 for r in rows:
  seq=contigs[r['target']];taln=re.sub('[^A-Za-z]','',r['taln']).upper();cands=[]
  for strand,frame,s,e,aa in six_segments(seq,45):
   if taln and taln in aa:cands.append((len(aa),strand,frame,s,e,aa))
  if not cands:
   # very short/gapped alignment can fail exact substring; use genomic alignment interval containment as fallback
   try:
    x1=int(float(r['tstart']));x2=int(float(r['tend']));lo=min(x1,x2)-1;hi=max(x1,x2)
    for strand,frame,s,e,aa in six_segments(seq,45):
     if s<=lo and e>=hi:cands.append((len(aa),strand,frame,s,e,aa))
   except:pass
  if not cands:continue
  cands.sort(key=lambda x:x[0]);L,strand,frame,s,e,aa=cands[0];hsh=hashlib.sha256(aa.encode()).hexdigest(); rid='EXT_'+hsh[:16]
  recovered[hsh]={'orf_id':rid,'aa_sha256':hsh,'aa_sequence':aa,'aa_length':L,'target_contig':r['target'],'source_accession':r['target'].split('|',1)[0],'strand':strand,'frame':frame,'start0':s,'end0':e}
  exact=(hsh==qt[r['query']]['aa_sha256'])
  qlen=int(qt[r['query']]['aa_length']);ratio=L/qlen
  core=(not exact and float(r['evalue'])<=1e-5 and float(r['qcov'])>=0.70 and float(r['fident'])>=0.25)
  broad=(not exact and float(r['evalue'])<=1e-3 and float(r['qcov'])>=0.50)
  longer=(L>=qlen+30 and ratio>=1.25)
  major=(longer and L>=max(180,1.5*qlen))
  annot.append({**r,'recovered_orf_id':rid,'recovered_aa_sha256':hsh,'recovered_aa_length':L,'source_accession':r['target'].split('|',1)[0],'orf_strand':strand,'orf_frame':frame,'orf_start0':s,'orf_end0':e,'exact_aa':exact,'core_homolog':core,'broad_homolog':broad,'longer_homolog':longer,'major_extension':major,'length_ratio':ratio,'component_id':comp[r['query']]})
 fa=out/'STAGE05F3A_RETURNED_RECOVERED_ORFS.faa'
 with fa.open('w') as f:
  for hsh,x in sorted(recovered.items()):f.write(f'>{x["orf_id"]} sha256={hsh} len={x["aa_length"]} source={x["source_accession"]}\n{x["aa_sequence"]}\n')
 dom=out/'STAGE05F3A_RETURNED_ORFS_RVT_ALL.domtblout';tbl=out/'STAGE05F3A_RETURNED_ORFS_RVT_ALL.tblout'
 subprocess.run(['hmmscan','--cpu',str(a.threads),'-E',str(EXPLORE),'--domE',str(EXPLORE),'--domtblout',str(dom),'--tblout',str(tbl),a.hmm,str(fa)],check=True,stdout=(out/'hmmscan.stdout').open('w'),stderr=(out/'hmmscan.stderr').open('w'))
 hb=defaultdict(list)
 for x in parse_domtbl(dom):hb[x['query']].append(x)
 byorf={x['orf_id']:x for x in recovered.values()}
 hmmrows=[]
 for oid,x in byorf.items():
  hits=sorted(hb.get(oid,[]),key=lambda z:(z['full_e'],z['i_e']));b=hits[0] if hits else None
  formal=any(z['full_e']<=FORMAL and z['i_e']<=FORMAL for z in hits);expl=bool(hits)
  hmmrows.append({'orf_id':oid,'aa_sha256':x['aa_sha256'],'aa_length':x['aa_length'],'source_accession':x['source_accession'],'formal_rvt_profile':formal,'exploratory_rvt_profile':expl,'best_rvt_target':b['hmm'] if b else '','best_full_evalue':b['full_e'] if b else '','best_domain_i_evalue':b['i_e'] if b else '','query_coverage':((b['ali_to']-b['ali_from']+1)/b['query_len']) if b else '','model_coverage':((b['hmm_to']-b['hmm_from']+1)/b['hmm_len']) if b else '','yxdd_like':bool(re.search(r'[YF][A-Z]DD',x['aa_sequence']))})
 hm={r['orf_id']:r for r in hmmrows}
 for r in annot:
  z=hm[r['recovered_orf_id']];r.update({k:z[k] for k in ['formal_rvt_profile','exploratory_rvt_profile','best_rvt_target','best_full_evalue','best_domain_i_evalue','query_coverage','model_coverage','yxdd_like']})
 af=out/'STAGE05F3A_ALL_HITS_ANNOTATED.tsv';cols=list(annot[0]) if annot else []
 with af.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(annot)
 hf=out/'STAGE05F3A_RECOVERED_ORF_HMM.tsv'
 with hf.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(hmmrows[0]) if hmmrows else ['orf_id'],delimiter='\t');w.writeheader();w.writerows(hmmrows)
 # summarize query and component endpoints
 qsum=[]
 for q in sorted(qseq):
  rr=[r for r in annot if r['query']==q];core=[r for r in rr if r['core_homolog']];long=[r for r in core if r['longer_homolog']];formal=[r for r in long if r['formal_rvt_profile']];maj=[r for r in formal if r['major_extension']]
  accs=sorted(set(r['source_accession'] for r in formal));best=min(formal,key=lambda r:float(r['best_full_evalue'])) if formal else None
  qsum.append({'query_id':q,'component_id':comp[q],'query_aa_length':qt[q]['aa_length'],'external_hit_n':len(rr),'core_homolog_n':len(core),'longer_core_homolog_n':len(long),'formal_longer_homolog_n':len(formal),'formal_major_extension_n':len(maj),'formal_longer_accession_n':len(accs),'formal_longer_accessions':';'.join(accs),'best_formal_orf':best['recovered_orf_id'] if best else '','best_formal_rvt_target':best['best_rvt_target'] if best else '','best_formal_evalue':best['best_full_evalue'] if best else '','best_formal_model_coverage':best['model_coverage'] if best else '','best_formal_query_coverage':best['query_coverage'] if best else '','best_formal_yxdd_like':best['yxdd_like'] if best else ''})
 qf=out/'STAGE05F3A_QUERY_SUMMARY.tsv'
 with qf.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(qsum[0]),delimiter='\t');w.writeheader();w.writerows(qsum)
 comps=defaultdict(list)
 for r in qsum:comps[r['component_id']].append(r)
 csum=[]
 for c,rr in sorted(comps.items()):
  csum.append({'component_id':c,'query_n':len(rr),'query_ids':';'.join(x['query_id'] for x in rr),'formal_longer_positive':any(int(x['formal_longer_homolog_n'])>0 for x in rr),'formal_major_positive':any(int(x['formal_major_extension_n'])>0 for x in rr),'formal_longer_accession_n':len(set(a for x in rr for a in x['formal_longer_accessions'].split(';') if a))})
 cf=out/'STAGE05F3A_COMPONENT_SUMMARY.tsv'
 with cf.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(csum[0]),delimiter='\t');w.writeheader();w.writerows(csum)
 formal_rows=[r for r in annot if r['core_homolog'] and r['longer_homolog'] and r['formal_rvt_profile']]
 ff=out/'STAGE05F3A_FORMAL_LONGER_HOMOLOGS.tsv'
 with ff.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(formal_rows)
 gate={'stage':'Stage05F3A1','gate':'PASS_F3A1_FROZEN_QUERY_EXTERNAL_EXPANSION_COMPLETE','gate_pass':True,'query_n':30,'component_n':26,'mmseqs_hit_rows_n':len(rows),'returned_recovered_unique_orf_n':len(recovered),'external_core_expansion_query_n':sum(int(x['core_homolog_n'])>0 for x in qsum),'external_longer_query_n':sum(int(x['longer_core_homolog_n'])>0 for x in qsum),'external_formal_longer_query_n':sum(int(x['formal_longer_homolog_n'])>0 for x in qsum),'external_formal_longer_component_n':sum(bool(x['formal_longer_positive']) for x in csum),'external_formal_major_extension_query_n':sum(int(x['formal_major_extension_n'])>0 for x in qsum),'external_formal_longer_multirun_query_n':sum(int(x['formal_longer_accession_n'])>=2 for x in qsum),'primary_endpoint':'external_formal_longer_component_n','claim_boundary':'External longer homolog plus frozen formal RVT support is RT-like evidence only; it does not establish biochemical RT activity or ART identity.'}
 (out/'STAGE05F3A_GATE.json').write_text(json.dumps(gate,indent=2)+'\n');print(json.dumps(gate,indent=2))
if __name__=='__main__':main()
