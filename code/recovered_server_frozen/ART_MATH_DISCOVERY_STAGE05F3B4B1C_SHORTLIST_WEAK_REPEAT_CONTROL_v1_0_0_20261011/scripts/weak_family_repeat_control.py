#!/usr/bin/env python3
import csv,json,hashlib,math,re,statistics,sys
from pathlib import Path
from collections import defaultdict,Counter

# Exact copy of B4B1 detector constants/logic.
RCMAP=str.maketrans('ACGTRYMKSWBDHVNacgtrymkswbdhvn','TGCAYRKMSWVHDBNtgcayrkm swvhd bn'.replace(' ',''))
SEED_K=(8,10,12,14,16)
MIN_COPIES=3
GAP_MIN=60
GAP_MAX=360
MAX_MOTIF=80
SCAN_FLANK=10000
ABSENCE_UPSTREAM_COMPLETE=5000
CODON={
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

def sha256_file(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(16<<20),b''): h.update(b)
 return h.hexdigest()
def clean_entry(b): return b.rstrip(b'\x00\r\n')
def revcomp(s): return s.translate(RCMAP)[::-1]
def translate(seq,frame): return ''.join(CODON.get(seq[i:i+3],'X') for i in range(frame,len(seq)-2,3))
def fetch_db(dbfh,off,ln):
 dbfh.seek(off); b=dbfh.read(ln)
 if len(b)!=ln: raise SystemExit(f'FAIL short DB read off={off} got={len(b)} expected={ln}')
 return re.sub(r'\s+','',clean_entry(b).decode('ascii','strict')).upper()
def load_index(path,wanted):
 out={}
 with open(path,encoding='ascii') as f:
  for line in f:
   if not line.strip(): continue
   x=line.split()
   if len(x)<3: raise SystemExit(f'FAIL malformed index {line[:120]!r}')
   if x[0] in wanted: out[x[0]]=(int(x[1]),int(x[2]))
 miss=wanted-set(out)
 if miss: raise SystemExit(f'FAIL missing DB keys n={len(miss)} examples={sorted(miss)[:10]}')
 return out
def exact_hits(seq,pep):
 hits=[]
 for f in (0,1,2):
  aa=translate(seq,f); pos=aa.find(pep)
  while pos>=0:
   s=f+3*pos; hits.append(('+',f,s,s+3*len(pep))); pos=aa.find(pep,pos+1)
 rc=revcomp(seq); n=len(seq)
 for f in (0,1,2):
  aa=translate(rc,f); pos=aa.find(pep)
  while pos>=0:
   rs=f+3*pos; re_=rs+3*len(pep); hits.append(('-',f,n-re_,n-rs)); pos=aa.find(pep,pos+1)
 return hits
def entropy(s):
 n=len(s)
 if not n:return 0.0
 c=Counter(s); return -sum((v/n)*math.log2(v/n) for v in c.values())
def low_complexity(s):
 if not s:return True
 return entropy(s)<1.35 or max(Counter(s).values())/len(s)>0.72
def chains_from_positions(pos):
 pos=sorted(pos); chains=[]; cur=[]
 for p in pos:
  if not cur: cur=[p]
  elif GAP_MIN<=p-cur[-1]<=GAP_MAX: cur.append(p)
  else:
   if len(cur)>=MIN_COPIES:chains.append(cur)
   cur=[p]
 if len(cur)>=MIN_COPIES:chains.append(cur)
 return chains
def extend_exact_core(seq,positions,k):
 left=0;right=k
 while left>-MAX_MOTIF and -left+right<MAX_MOTIF:
  idx=[p+left-1 for p in positions]
  if min(idx)<0:break
  if len({seq[i] for i in idx})!=1:break
  left-=1
 while -left+right<MAX_MOTIF:
  idx=[p+right for p in positions]
  if max(idx)>=len(seq):break
  if len({seq[i] for i in idx})!=1:break
  right+=1
 starts=[p+left for p in positions]
 return seq[starts[0]:starts[0]+(-left+right)],starts
def scan_region(seq,region_label,offset_rel_candidate):
 raw=[]
 for k in SEED_K:
  d=defaultdict(list)
  for i in range(0,len(seq)-k+1):
   mer=seq[i:i+k]
   if 'N' in mer or low_complexity(mer):continue
   d[mer].append(i)
  for mer,pos in d.items():
   if len(pos)<MIN_COPIES:continue
   for chain in chains_from_positions(pos):
    motif,starts=extend_exact_core(seq,chain,k)
    if len(motif)<8 or low_complexity(motif):continue
    gaps=[starts[i+1]-starts[i] for i in range(len(starts)-1)]
    spacers=[g-len(motif) for g in gaps]
    if any(x<0 for x in spacers):continue
    mean=statistics.mean(spacers) if spacers else 0.0
    med=statistics.median(spacers) if spacers else 0.0
    sd=statistics.pstdev(spacers) if len(spacers)>1 else 0.0
    cv=(sd/mean) if mean>0 else 0.0
    rel=[offset_rel_candidate+x for x in starts]
    strict=(region_label=='UPSTREAM' and 3<=len(starts)<=21 and all(120<=x<=220 for x in spacers))
    broad=(region_label=='UPSTREAM' and 3<=len(starts)<=30 and 100<=med<=240 and (sum(80<=x<=280 for x in spacers)/len(spacers) if spacers else 0)>=0.75)
    raw.append({'region':region_label,'seed_k':k,'motif':motif,'motif_len':len(motif),'copy_n':len(starts),
      'local_repeat_starts':','.join(map(str,starts)),'candidate_relative_repeat_starts':','.join(map(str,rel)),
      'start_to_start_gaps':','.join(map(str,gaps)),'spacer_lengths':','.join(map(str,spacers)),
      'median_spacer_bp':med,'min_spacer_bp':min(spacers) if spacers else '', 'max_spacer_bp':max(spacers) if spacers else '',
      'spacer_cv':cv,'array_span_bp':starts[-1]+len(motif)-starts[0],
      'strict_published_art_geometry':strict,'broad_art_like_geometry':broad})
 raw.sort(key=lambda r:(-r['copy_n'],-r['motif_len'],r['spacer_cv'],r['motif']))
 keep=[]
 for r in raw:
  s=set(map(int,r['local_repeat_starts'].split(','))); dup=False
  for q in keep:
   t=set(map(int,q['local_repeat_starts'].split(','))); denom=min(len(s),len(t))
   if denom and len(s&t)/denom>=0.8 and abs(r['motif_len']-q['motif_len'])<=12:dup=True;break
  if not dup:keep.append(r)
 return keep

def main():
 if len(sys.argv)!=12:
  raise SystemExit('usage: weak_family_repeat_control.py DSUM FM IM ORFDB ORFIDX GDB GIDX B4B1_FAMSUM OUTDIR LOCK N_EXPECT')
 dsump,fmp,imp,orfdb,orfidx,gdb,gidx,b4b1sum,outdir,lockp,nexp=sys.argv[1:]
 out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
 lock=json.load(open(lockp));assert lock['gate_pass'] is True
 drows=list(csv.DictReader(open(dsump),delimiter='\t'))
 weak=[r for r in drows if r['final_best_tier']=='WEAK']
 unresolved=[r for r in drows if r['final_best_tier'].startswith('BOLTZ_TECHNICAL_FAILURE')]
 if len(drows)!=32 or len(weak)!=9 or len(unresolved)!=7: raise SystemExit(f'FAIL DSUM 32={len(drows)} weak={len(weak)} unresolved={len(unresolved)}')
 weak_fams=[r['family_rep_id'] for r in weak]; weakset=set(weak_fams)
 tier={r['family_rep_id']:r['final_best_tier'] for r in weak}
 pri={r['family_rep_id']:r['candidate_id'] for r in weak}
 members=[]
 with open(fmp) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['family_rep_id'] in weakset:members.append(r)
 im={}
 wanted={r['candidate_id'] for r in members}
 with open(imp) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['instance_id'] in wanted:im[r['instance_id']]=r
 if wanted-set(im):raise SystemExit('FAIL missing provenance')
 orfkeys={im[r['candidate_id']]['orfdb_key'] for r in members};contigkeys={im[r['candidate_id']]['parent_contig_key'] for r in members}
 oi=load_index(orfidx,orfkeys);gi=load_index(gidx,contigkeys)
 peps={};seqs={}
 with open(orfdb,'rb',buffering=0) as fh:
  for k,(o,l) in oi.items():peps[k]=fetch_db(fh,o,l).rstrip('*')
 with open(gdb,'rb',buffering=0) as fh:
  for k,(o,l) in gi.items():seqs[k]=fetch_db(fh,o,l)
 arrays=[];locus=[]
 for ix,m in enumerate(members,1):
  q=m['candidate_id'];p=im[q];pep=peps[p['orfdb_key']];seq=seqs[p['parent_contig_key']]
  hits=exact_hits(seq,pep)
  if len(hits)!=1: raise SystemExit(f'FAIL control exact placement {q} n={len(hits)}')
  strand,frame,s,e=hits[0]
  if strand=='+':
   up=seq[max(0,s-SCAN_FLANK):s]; dn=seq[e:min(len(seq),e+SCAN_FLANK)]
  else:
   up=revcomp(seq[e:min(len(seq),e+SCAN_FLANK)]); dn=revcomp(seq[max(0,s-SCAN_FLANK):s])
  up_av=len(up);dn_av=len(dn)
  up_arr=scan_region(up,'UPSTREAM',-up_av);dn_arr=scan_region(dn,'DOWNSTREAM',len(pep)*3);arr=up_arr+dn_arr
  for j,a in enumerate(arr,1):
   z={'family_rep_id':m['family_rep_id'],'structural_tier':'WEAK','candidate_id':q,'accession':p['rescued_accession'],'parent_contig_key':p['parent_contig_key'],'strand':strand,'upstream_available_bp':up_av,'downstream_available_bp':dn_av,'array_rank_within_locus':j};z.update(a);arrays.append(z)
  locus.append({'family_rep_id':m['family_rep_id'],'structural_tier':'WEAK','candidate_id':q,'accession':p['rescued_accession'],'upstream_available_bp':up_av,'downstream_available_bp':dn_av,
    'any_generic_repeat_array':bool(arr),'upstream_repeat_array_n':len(up_arr),'strict_published_art_array_n':sum(bool(a['strict_published_art_geometry']) for a in arr),'broad_art_like_array_n':sum(bool(a['broad_art_like_geometry']) for a in arr)})
  if ix%100==0:print(f'[SCAN] {ix}/{len(members)} arrays={len(arrays)}',flush=True)
 if int(nexp)>0 and len(members)!=int(nexp): print(f'[WARN] expected_n_hint={nexp} observed={len(members)}')
 # outputs
 arrp=out/'F3B4B1C_WEAK_ALL_DETECTED_REPEAT_ARRAYS.tsv';locp=out/'F3B4B1C_WEAK_LOCUS_SUMMARY.tsv'
 ac=list(arrays[0]) if arrays else ['family_rep_id','structural_tier','candidate_id','accession','parent_contig_key','strand','upstream_available_bp','downstream_available_bp','array_rank_within_locus','region','seed_k','motif','motif_len','copy_n','local_repeat_starts','candidate_relative_repeat_starts','start_to_start_gaps','spacer_lengths','median_spacer_bp','min_spacer_bp','max_spacer_bp','spacer_cv','array_span_bp','strict_published_art_geometry','broad_art_like_geometry']
 with arrp.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=ac,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(arrays)
 lc=list(locus[0]);
 with locp.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=lc,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(locus)
 fs=[]
 for fam in weak_fams:
  lr=[x for x in locus if x['family_rep_id']==fam]; ar=[x for x in arrays if x['family_rep_id']==fam]
  strict=[x for x in ar if x['strict_published_art_geometry']];broad=[x for x in ar if x['broad_art_like_geometry']];up=[x for x in ar if x['region']=='UPSTREAM']
  fs.append({'family_rep_id':fam,'structural_tier':'WEAK','priority_candidate_id':pri[fam],'family_instance_n':len(lr),'family_accession_n':len({x['accession'] for x in lr}),
    'instance_n_with_any_repeat_array':sum(x['any_generic_repeat_array'] for x in lr),'accession_n_with_any_repeat_array':len({x['accession'] for x in lr if x['any_generic_repeat_array']}),
    'instance_n_with_upstream_repeat_array':len({x['candidate_id'] for x in up}),'accession_n_with_upstream_repeat_array':len({x['accession'] for x in up}),
    'strict_published_art_array_n':len(strict),'strict_published_art_instance_n':len({x['candidate_id'] for x in strict}),'strict_published_art_accession_n':len({x['accession'] for x in strict}),
    'broad_art_like_array_n':len(broad),'broad_art_like_instance_n':len({x['candidate_id'] for x in broad}),'broad_art_like_accession_n':len({x['accession'] for x in broad})})
 fsp=out/'F3B4B1C_WEAK_FAMILY_REPEAT_SUMMARY_9.tsv';fc=list(fs[0])
 with fsp.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fc,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(fs)
 # survivor B4B1 summary
 sv=list(csv.DictReader(open(b4b1sum),delimiter='\t'))
 assert len(sv)==16
 comp=[]
 for metric in ['instance_n_with_any_repeat_array','instance_n_with_upstream_repeat_array','strict_published_art_instance_n','broad_art_like_instance_n']:
  svnum=sum(int(x[metric]) for x in sv);svden=sum(int(x['family_instance_n']) for x in sv)
  wknum=sum(int(x[metric]) for x in fs);wkden=sum(int(x['family_instance_n']) for x in fs)
  comp.append({'metric':metric,'survivor_positive_n':svnum,'survivor_denominator_n':svden,'survivor_rate':svnum/svden if svden else 0,
    'weak_positive_n':wknum,'weak_denominator_n':wkden,'weak_rate':wknum/wkden if wkden else 0,'rate_ratio_survivor_over_weak':(svnum/svden)/(wknum/wkden) if svden and wkden and wknum else 'INF' if svnum else 0})
 for metric in ['instance_n_with_any_repeat_array','instance_n_with_upstream_repeat_array','strict_published_art_array_n','broad_art_like_array_n']:
  svf=sum(int(x[metric])>0 for x in sv);wkf=sum(int(x[metric])>0 for x in fs)
  comp.append({'metric':'family_any__'+metric,'survivor_positive_n':svf,'survivor_denominator_n':16,'survivor_rate':svf/16,
    'weak_positive_n':wkf,'weak_denominator_n':9,'weak_rate':wkf/9,'rate_ratio_survivor_over_weak':(svf/16)/(wkf/9) if wkf else 'INF' if svf else 0})
 cp=out/'F3B4B1C_SURVIVOR_VS_WEAK_DESCRIPTIVE_COMPARISON.tsv';cc=list(comp[0])
 with cp.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=cc,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(comp)
 gate={'stage':'Stage05F3B4B1C','gate':'PASS_F3B4B1C_SAME_SHORTLIST_WEAK_REPEAT_CONTROL_COMPLETE','gate_pass':True,
  'role':'post-discovery matched structural comparator; same frozen 32-family priority panel','survivor_family_n':16,'weak_family_n':9,'technical_unresolved_family_n':7,
  'weak_member_locus_n':len(locus),'weak_detected_array_n':len(arrays),'weak_family_n_with_any_repeat_array':sum(x['instance_n_with_any_repeat_array']>0 for x in fs),
  'weak_family_n_with_upstream_repeat_array':sum(x['instance_n_with_upstream_repeat_array']>0 for x in fs),'weak_family_n_with_strict_geometry':sum(x['strict_published_art_array_n']>0 for x in fs),'weak_family_n_with_broad_geometry':sum(x['broad_art_like_array_n']>0 for x in fs),
  'control_definition_frozen_before_control_results':True,'control_was_prespecified_before_b4b1_survivor_results':False,
  'detector_parameters_changed':False,'unresolved_treated_as_negative':False,'novelty_decision_performed':False,
  'weak_family_summary_sha256':sha256_file(fsp),'comparison_sha256':sha256_file(cp),
  'claim_boundary':'Descriptive post-discovery same-shortlist comparator only. WEAK structural compatibility is not biochemical non-RT; unresolved families are censored, not negative.'}
 gp=out/'STAGE05F3B4B1C_GATE.json';gp.write_text(json.dumps(gate,indent=2)+'\n')
 print(json.dumps(gate,indent=2));print('\n=== SURVIVOR VS WEAK DESCRIPTIVE COMPARISON ===');print('\t'.join(cc));
 for x in comp:print('\t'.join(str(x[c]) for c in cc))
 print('\n=== WEAK FAMILY SUMMARY 9 ===');print('\t'.join(fc));
 for x in fs:print('\t'.join(str(x[c]) for c in fc))
 print('PASS_STAGE05F3B4B1C_SAME_SHORTLIST_WEAK_REPEAT_CONTROL_COMPLETE')
if __name__=='__main__':main()
