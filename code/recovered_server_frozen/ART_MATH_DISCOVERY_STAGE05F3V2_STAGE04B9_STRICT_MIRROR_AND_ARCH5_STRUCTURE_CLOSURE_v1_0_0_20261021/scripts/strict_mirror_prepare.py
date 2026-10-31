#!/usr/bin/env python3
import argparse,csv,hashlib,json,re,sys
from collections import defaultdict,Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from rule_core import repkey,select_families,qualifying_profile_pass

MOTIF_RE=re.compile(r'(?=([YF][A-Z]DD))')
EXPECTED_FORMAL=88487
EXPECTED_UNIQUE=85503
EXPECTED_FAMILIES=26396
EXPECTED_ACCESSIONS=254
EXPECTED_HIST_SCRIPT_SHA='d222f638c627566909765012874c7aaa531a1a79b5e67a0d73f726f184a30ed8'

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def shaseq(s):return hashlib.sha256(s.encode('ascii')).hexdigest()
def read_tsv(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,cols):
 p.parent.mkdir(parents=True,exist_ok=True)
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def fasta_iter(p):
 h=None;seq=[]
 with open(p,encoding='ascii') as f:
  for line in f:
   line=line.strip()
   if not line:continue
   if line.startswith('>'):
    if h is not None:yield h,''.join(seq).rstrip('*').upper()
    h=line[1:].split()[0];seq=[]
   else:seq.append(line)
  if h is not None:yield h,''.join(seq).rstrip('*').upper()
def one_dir(P,token):
 ds=sorted(d for d in P.iterdir() if d.is_dir() and token in d.name)
 if not ds:raise RuntimeError('cannot locate '+token)
 return ds[-1]

def load_native_delta_partition(D):
 em=D/'outputs/F3V2_FORMAL_EXACT_MAP.tsv'
 cl=D/'outputs/F3V2_FAMILY_cluster.tsv'
 gate=D/'outputs/STAGE05F3V2_GATE.json'
 for p in (em,cl,gate):
  if not p.is_file():raise RuntimeError('missing native Stage03 partition input '+str(p))
 g=json.load(open(gate))
 if g.get('gate_pass') is not True:raise RuntimeError('delta gate not PASS')
 if int(g.get('formal_rvt_instance_n',-1))!=EXPECTED_FORMAL or int(g.get('exact_aa_unique_n',-1))!=EXPECTED_UNIQUE or int(g.get('operational_family_n',-1))!=EXPECTED_FAMILIES:
  raise RuntimeError('delta Stage03 counts differ from frozen 88487/85503/26396')
 if g.get('exact_map_sha256') and sha(em)!=g['exact_map_sha256']:
  raise RuntimeError('native exact map SHA differs from Stage03 gate')
 exact_rows=read_tsv(em)
 if len(exact_rows)!=EXPECTED_FORMAL:raise RuntimeError(f'exact map rows={len(exact_rows)}')
 q2u={r['formal_query_id']:r['exact_unique_id'] for r in exact_rows}
 u2sha={r['exact_unique_id']:r['aa_sha256'] for r in exact_rows}
 if len(set(q2u.values()))!=EXPECTED_UNIQUE:raise RuntimeError('exact unique count mismatch')
 member_to_rep={}
 reps=set()
 with open(cl) as f:
  for line in f:
   if not line.strip():continue
   z=line.rstrip('\n').split('\t')
   if len(z)<2:continue
   rep,mem=z[0],z[1]
   if mem in member_to_rep and member_to_rep[mem]!=rep:raise RuntimeError('native exact member maps to >1 family')
   member_to_rep[mem]=rep;reps.add(rep)
 allu=set(q2u.values())
 missing=allu-set(member_to_rep)
 if missing:raise RuntimeError(f'native cluster membership missing exact IDs n={len(missing)} examples={sorted(missing)[:10]}')
 if len(reps)!=EXPECTED_FAMILIES:raise RuntimeError(f'native family n={len(reps)} expected={EXPECTED_FAMILIES}')
 q2fam={q:member_to_rep[u] for q,u in q2u.items()}
 return q2u,u2sha,q2fam,{'exact_map':str(em),'cluster_tsv':str(cl),'gate':str(gate),'exact_map_sha256':sha(em),'cluster_tsv_sha256':sha(cl),'gate_sha256':sha(gate)}

def discover_existing_deltaonly(B2):
 t=B2/'outputs/F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.tsv'
 f=B2/'outputs/F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.faa'
 if not t.is_file() or not f.is_file():return {},{}
 rows=read_tsv(t); seq=dict(fasta_iter(f))
 bysha={shaseq(s):cid for cid,s in seq.items()}
 return {r['priority_candidate_id']:r for r in rows},bysha

def discover_b3_cif_by_sequence(B3,B2):
 # Exact-AA reuse only. Candidate identity or exact sequence equality is required.
 fa=B2/'outputs/F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.faa'
 if not fa.is_file():return {}
 seq=dict(fasta_iter(fa)); cid2sha={cid:shaseq(s) for cid,s in seq.items()}
 # Read any status TSV carrying candidate_id/cif_path; prefer paths that still exist.
 cidpaths=defaultdict(list)
 for p in B3.rglob('*.tsv'):
  try:
   with open(p,newline='') as f:
    r=csv.DictReader(f,delimiter='\t'); fields=set(r.fieldnames or [])
    if not {'candidate_id','cif_path'}.issubset(fields):continue
    for row in r:
     cid=row.get('candidate_id',''); cp=row.get('cif_path','')
     if cid in cid2sha and cp and Path(cp).is_file():cidpaths[cid].append(Path(cp))
  except Exception:pass
 # Fallback direct recursive candidate directory search.
 for cid in cid2sha:
  if cidpaths[cid]:continue
  for p in B3.rglob('*_model_0.cif'):
   if cid in str(p):cidpaths[cid].append(p)
 out={}
 bysha=defaultdict(list)
 for cid,paths in cidpaths.items():
  uniq={sha(p):p for p in paths if p.is_file()}
  if len(uniq)==1:bysha[cid2sha[cid]].append(next(iter(uniq.values())))
 for s,paths in bysha.items():
  # If several candidate IDs with same AA have structures, require all file hashes identical; otherwise rerun.
  uniq={sha(p):p for p in paths}
  if len(uniq)==1:out[s]=next(iter(uniq.values()))
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--project-root',required=True);ap.add_argument('--home',required=True);ap.add_argument('--semantics-audit',required=True)
 a=ap.parse_args();P=Path(a.project_root).resolve();H=Path(a.home).resolve();O=H/'outputs';I=H/'inputs';W=H/'work';O.mkdir(exist_ok=True);I.mkdir(exist_ok=True);W.mkdir(exist_ok=True)
 D=one_dir(P,'STAGE05F3V2_60G_DELTA_DISTRIBUTED_ACQUISITION_R2_v1_0_2_20261013')
 B2=one_dir(P,'STAGE04B2_HISTORICAL_F3B3A_EXACT_PRIORITY_PROJECTION')
 B3=one_dir(P,'STAGE04B3_DELTA32_BOLTZ2_FROZEN_R5C')
 sa=json.load(open(a.semantics_audit));
 if sa.get('gate_pass') is not True:raise RuntimeError('historical semantics audit not PASS')
 qdef=sa.get('qualifying_definition_used_for_delta');
 if not qdef:raise RuntimeError('missing qualifying definition')
 q2u,u2sha,q2fam,partprov=load_native_delta_partition(D)
 finst=D/'outputs/F3V2_FORMAL_INSTANCES.tsv'; manp=D/'frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'
 formal_byacc=defaultdict(dict);formal_n=0
 with open(finst,newline='') as f:
  r=csv.DictReader(f,delimiter='\t')
  for row in r:
   if row.get('formal')!='True':continue
   formal_byacc[row['source_accession']][row['query_id']]=row;formal_n+=1
 if formal_n!=EXPECTED_FORMAL:raise RuntimeError(f'formal n {formal_n}')
 mans=read_tsv(manp);accs=[r['accession'] for r in mans]
 if len(accs)!=EXPECTED_ACCESSIONS or len(set(accs))!=EXPECTED_ACCESSIONS:raise RuntimeError('accession universe changed')
 formal_acc=[x for x in accs if x in formal_byacc];zero_acc=[x for x in accs if x not in formal_byacc]
 if len(formal_acc)+len(zero_acc)!=EXPECTED_ACCESSIONS:raise RuntimeError('accession accounting')
 cand=[];byfam=defaultdict(list)
 for ai,acc in enumerate(formal_acc,1):
  rows=formal_byacc[acc]
  faa=D/'outputs/formal_faa'/f'{acc}.formal.faa';dom=D/'work/per_accession'/acc/'RVT_ALL.domtblout'
  if not faa.is_file() or not dom.is_file():raise RuntimeError('missing per-accession frozen output '+acc)
  seq={q:s for q,s in fasta_iter(faa)}
  if set(seq)!=set(rows):raise RuntimeError(f'formal FASTA mismatch {acc}')
  formal_prof=defaultdict(set);qual_prof=defaultdict(set);bestdom={}
  with open(dom,errors='ignore') as f:
   for line in f:
    if not line.strip() or line.startswith('#'):continue
    x=line.split();
    if len(x)<22:continue
    profile=x[0];q=x[3]
    if q not in rows:continue
    full=float(x[6]);ie=float(x[12]);ali_from=int(x[17]);ali_to=int(x[18]);env_from=int(x[19]);env_to=int(x[20])
    if full<=1e-4 and ie<=1e-4:
     formal_prof[q].add(profile);rec={'profile':profile,'full_evalue':full,'domain_i_evalue':ie,'ali_from':ali_from,'ali_to':ali_to,'env_from':env_from,'env_to':env_to}
     key=(ie,full,profile,ali_from,ali_to);old=bestdom.get(q)
     if old is None or key < (old['domain_i_evalue'],old['full_evalue'],old['profile'],old['ali_from'],old['ali_to']):bestdom[q]=rec
    if qualifying_profile_pass(qdef,full,ie):qual_prof[q].add(profile)
  for q,row in rows.items():
   s=seq[q];u=q2u.get(q);fr=q2fam.get(q)
   if u is None or fr is None:raise RuntimeError('native partition missing formal candidate '+q)
   d=bestdom.get(q)
   if d is None:raise RuntimeError('formal candidate lacks formal domain '+q)
   if float(row['best_domain_i_evalue'])!=float(d['domain_i_evalue']):raise RuntimeError('best-domain mismatch '+q)
   lo=max(1,d['env_from']);hi=min(len(s),d['env_to']);sub=s[lo-1:hi]
   hits=[(lo+m.start(),m.group(1)) for m in MOTIF_RE.finditer(sub)];gh=[(1+m.start(),m.group(1)) for m in MOTIF_RE.finditer(s)]
   c={'candidate_id':q,'original_query_name':q,'family_rep_id':fr,'exact_rep_id':u,'accession':acc,'aa_length':len(s),'formal_model_n':len(formal_prof[q]),'qualifying_model_n':len(qual_prof[q]),'best_full_evalue':float(row['full_evalue']),'best_domain_i_evalue':float(row['best_domain_i_evalue']),'best_formal_profile':d['profile'],'best_domain_env_from':lo,'best_domain_env_to':hi,'motif_c_like_in_best_domain':bool(hits),'motif_c_like_domain_n':len(hits),'motif_c_like_domain_first_pos':hits[0][0] if hits else '','motif_c_like_domain_first_seq':hits[0][1] if hits else '','motif_c_like_anywhere':bool(gh),'motif_c_like_anywhere_n':len(gh),'sequence_sha256':shaseq(s)}
   cand.append(c);byfam[fr].append(c)
  if ai%25==0 or ai==len(formal_acc):print(f'[STRICT MIRROR] accession {ai}/{len(formal_acc)} candidates={len(cand)}',flush=True)
 if len(cand)!=EXPECTED_FORMAL or len(byfam)!=EXPECTED_FAMILIES:raise RuntimeError(f'mirror coverage candidates={len(cand)} fam={len(byfam)}')
 famrows=[]
 for fr in sorted(byfam):
  m=byfam[fr];rep=sorted(m,key=repkey)[0];motif_n=sum(bool(x['motif_c_like_in_best_domain']) for x in m)
  famrows.append({'family_rep_id':fr,'candidate_instance_n':len(m),'accession_n':len({x['accession'] for x in m}),'multirun_replication':len({x['accession'] for x in m})>=2,'max_aa_length':max(x['aa_length'] for x in m),'max_formal_model_n':max(x['formal_model_n'] for x in m),'max_qualifying_model_n':max(x['qualifying_model_n'] for x in m),'motif_c_like_candidate_n':motif_n,'motif_c_like_family':motif_n>0,'priority_candidate_id':rep['candidate_id'],'priority_candidate_accession':rep['accession'],'priority_candidate_aa_length':rep['aa_length'],'priority_candidate_formal_model_n':rep['formal_model_n'],'priority_candidate_best_domain_i_evalue':rep['best_domain_i_evalue'],'priority_candidate_motif_c_like':rep['motif_c_like_in_best_domain'],'priority_candidate_best_formal_profile':rep['best_formal_profile']})
 selected,reason,rank,fam_by_id=select_families(famrows)
 for x in famrows:
  fr=x['family_rep_id'];x['shortlisted']=fr in rank;x['shortlist_rank']=rank.get(fr,'');x['shortlist_reason']=';'.join(reason.get(fr,[]))
 famrows.sort(key=lambda x:(0 if x['shortlisted'] else 1,x['shortlist_rank'] if x['shortlisted'] else 10**9,-x['accession_n'],-x['max_aa_length'],x['family_rep_id']))
 write_tsv(O/'STRICT_MIRROR_CANDIDATE_CHARACTERIZATION_88487.tsv',cand,list(cand[0]))
 write_tsv(O/'STRICT_MIRROR_FAMILY_PRIORITY_26396.tsv',famrows,list(famrows[0]))
 short=[];need=defaultdict(set)
 for fr in selected:
  x=dict(fam_by_id[fr]);x.update(shortlisted=True,shortlist_rank=rank[fr],shortlist_reason=';'.join(reason[fr]));short.append(x);need[x['priority_candidate_accession']].add(x['priority_candidate_id'])
 write_tsv(O/'STRICT_MIRROR_SHORTLIST_32.tsv',short,list(short[0]))
 seqsel={}
 for acc,qs in need.items():
  for q,s in fasta_iter(D/'outputs/formal_faa'/f'{acc}.formal.faa'):
   if q in qs:seqsel[q]=s
 if set(seqsel)!={x['priority_candidate_id'] for x in short}:raise RuntimeError('strict mirror FASTA recovery mismatch')
 with open(O/'STRICT_MIRROR_SHORTLIST_32.faa','w') as f:
  for x in short:
   q=x['priority_candidate_id'];s=seqsel[q];f.write(f'>{q} family={x["family_rep_id"]} accession={x["priority_candidate_accession"]} aa_length={len(s)}\n')
   for i in range(0,len(s),80):f.write(s[i:i+80]+'\n')
 # Compare to existing novelty-only DELTA32 and prepare strict-mirror structure input/reuse.
 _,b2sha=discover_existing_deltaonly(B2);b3reuse=discover_b3_cif_by_sequence(B3,B2)
 overlap=[];manifest=[];yd=I/'strict_mirror_boltz_yaml';yd.mkdir(parents=True,exist_ok=True);sd=O/'strict_mirror_structures';sd.mkdir(parents=True,exist_ok=True)
 for i,x in enumerate(short,1):
  qid=f'SMQ{i:04d}';cid=x['priority_candidate_id'];s=seqsel[cid];ss=shaseq(s);oldcid=b2sha.get(ss,'');reuse=b3reuse.get(ss)
  status='NEEDS_BOLTZ';reuse_path=''
  if reuse and reuse.is_file():
   import shutil;shutil.copy2(reuse,sd/f'{qid}.cif');status='REUSED_EXISTING_B3_EXACT_AA';reuse_path=str(reuse)
  y=yd/f'{qid}.yaml';y.write_text('version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: '+s+'\n')
  manifest.append({'query_id':qid,'branch':'STRICT_MIRROR','shortlist_rank':i,'family_rep_id':x['family_rep_id'],'candidate_id':cid,'accession':x['priority_candidate_accession'],'aa_length':len(s),'sequence_sha256':ss,'deltaonly32_exact_aa_match_candidate_id':oldcid,'structure_prepare_status':status,'reused_cif_source':reuse_path,'yaml_path':str(y.resolve())})
  overlap.append({'strict_mirror_rank':i,'strict_mirror_family_rep_id':x['family_rep_id'],'strict_mirror_candidate_id':cid,'sequence_sha256':ss,'deltaonly32_exact_aa_overlap':bool(oldcid),'deltaonly32_candidate_id':oldcid})
 write_tsv(O/'STRICT_MIRROR_STRUCTURE_MANIFEST_32.tsv',manifest,list(manifest[0]))
 write_tsv(O/'STRICT_MIRROR_VS_DELTAONLY32_EXACT_AA_OVERLAP.tsv',overlap,list(overlap[0]))
 gate={'stage':'Stage05F3V2-Stage04B9A','gate':'PASS_STAGE04B9A_STRICT_MIRROR_26396_TO_32_PREPARED','gate_pass':True,'delta_accession_n':EXPECTED_ACCESSIONS,'formal_bearing_accession_n':len(formal_acc),'zero_formal_accession_n':len(zero_acc),'formal_instance_n':len(cand),'exact_aa_unique_n':EXPECTED_UNIQUE,'native_operational_family_n':len(famrows),'family_universe':'ALL_NATIVE_DELTA_OPERATIONAL_FAMILIES_NO_COMBINED_NOVELTY_FILTER','historical_rule':'exact F3B3A 4x8 quotas plus deterministic backfill','shortlist_n':32,'deltaonly32_exact_aa_overlap_n':sum(bool(x['deltaonly32_exact_aa_overlap']) for x in overlap),'strict_mirror_structure_reused_n':sum(x['structure_prepare_status'].startswith('REUSED') for x in manifest),'strict_mirror_structure_missing_n':sum(x['structure_prepare_status']=='NEEDS_BOLTZ' for x in manifest),'native_partition_provenance':partprov,'qualifying_model_definition':qdef,'known_annotation_used_for_selection':False,'combined_old_delta_family_filter_used':False,'structure_used_for_selection':False,'retuning_performed':False,'science_change':False}
 with open(O/'STAGE04B9A_PREPARE_GATE.json','w') as f:json.dump(gate,f,indent=2);f.write('\n')
 print(json.dumps(gate,indent=2));print('PASS_STAGE04B9A_STRICT_MIRROR_26396_TO_32_PREPARED')
if __name__=='__main__':
 try:main()
 except Exception as e:
  print('[HOLD]',e,file=sys.stderr);raise
