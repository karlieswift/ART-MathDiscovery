#!/usr/bin/env python3
"""J2 external-cohort admission, prefreeze, scoring and postfreeze evaluation.

This coordinator uses the Python standard library. Actual scoring delegates to
byte-identical vendored project engines, which require numpy/pandas and HMMER.
It never manufactures scientific labels or substitutes a new scoring algorithm.
"""
from __future__ import annotations
import argparse, csv, datetime as dt, hashlib, importlib, importlib.metadata, json, math
import os, platform, random, re, resource, shutil, signal, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
METHODS = ['ARM_M_FULL_FROZEN', 'PRODIGAL_RVT_PLUS_FROZEN_REC',
           'SIXFRAME_RVT_PLUS_FROZEN_REC', 'PRODIGAL_RVT_ONLY', 'SIXFRAME_RVT_ONLY']
PRIMARY = METHODS[:3]
DEPENDENCIES = ['numpy','pandas','tqdm']
RVT_SHA = '74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e'
CANDIDATE_FIELDS = ['sample_id','accession','accession_aliases','native_parent_accession','genome_path',
 'fasta_sha256','sequence_sha256','taxon','source_url','source_version',
 'source_study','source_publication_date','family_cluster','leakage_unit',
 'match_group','exposure_status','eligibility_status','eligibility_reason',
 'leakage_review_path','leakage_review_sha256','leakage_review_status']
TRUTH_FIELDS = ['sample_id','target_id','label_class','label_scope','truth_start0','truth_end0',
 'array_start0','array_end0','rt_start0','rt_end0','author_truth_url',
 'truth_evidence_path','truth_evidence_sha256','truth_basis','curator','truth_fixed_utc']
POLICY = {
 'schema_version':'J2-1.0.0', 'analysis_kind':'independent_external_benchmark',
 'primary_methods':PRIMARY, 'secondary_methods':METHODS[3:],
 'minimum_sequence_length_bp':4096, 'allowed_alphabet':'ACGT',
 'multi_contig_policy':'one complete native contig per row; no concatenation; no truth-centred cropping',
 'windows_bp':[128,256,512,1024,2048,4096], 'stride_fraction':0.25,
 'strands':['+','-'], 'proposal_q_threshold':0.005,
 'arm_m_engine':'exact formal ACNV+recurrence orchestration of unchanged original helpers, followed by frozen clean_two_arm build_two_arm_board',
 'execution_mode':'exact_two_arm',
 'unused_legacy_irregularity':'omitted in exact_two_arm; available with full_legacy; never an input to formal two-arm board',
 'arm_m_ranking':['support_family_n desc','sum_two_family_surprisal desc','best_single_family_q asc','start0 asc'],
 'rvt_expected_sha256':RVT_SHA, 'rvt_full_evalue':1e-4,'rvt_domain_ievalue':1e-4,
 'rt_radius_bp':10000,'sixframe_min_aa':45,'prodigal_mode':'meta',
 'rt_locus_merge':'same route and strand, strict overlap, gap=0',
 'rt_plus_rec_ranking':['has_recurrence desc','best_local_recurrence_q asc','full_evalue asc','domain_i_evalue asc','hmm_score desc'],
 'rt_only_ranking':['full_evalue asc','domain_i_evalue asc','hmm_score desc'],
 'candidate_limit':None, 'inputs':'identical full contig for all methods',
 'primary_target':'independently author-defined ART system interval; strict genomic overlap, not functional validation',
 'component_targets':'report array localization and RT localization separately where independently annotated',
 'failure_policy':'all admitted samples retained; failed methods contribute zero recovery and remain flagged; burden undefined for failed methods',
 'negative_policy':'candidate burden only; background and hard-neighbor contigs are not certified whole-genome ART negatives',
 'exposure_policy':'report nominated-span fraction and RT-neighborhood fraction; never impose a truth-derived crop',
 'bootstrap_unit':'prefrozen leakage_unit; paired resampling across methods',
 'bootstrap_replicates':2000,'bootstrap_seed':20261022,'bootstrap_minimum_positive_units':5,
 'small_sample_policy':'below 5 independent positive leakage units: descriptive metrics only, no interval',
 'truth_used_for_scoring':False, 'retuning_authorized':False,
 'computational_blinding_scope':'script-level separation of truth from score inputs; no assertion of human/author blinding',
 'legacy_stage_labels':'unchanged source engine outputs retain historical Stage03A labels; cohort eligibility is defined only by this driver analysis_kind and admission evidence'
}

class Hold(Exception): pass
def now(): return dt.datetime.now(dt.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for x in iter(lambda:f.read(1024*1024),b''): h.update(x)
 return h.hexdigest()
def digest(x): return hashlib.sha256(x.encode()).hexdigest()
def read_tsv(p):
 with Path(p).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,fields=None):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 fields=fields or list(rows[0] if rows else {'empty':''})
 with p.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def write_json(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def js(p):return json.loads(Path(p).read_text())
def resolved(s,base): return (Path(s) if Path(s).is_absolute() else Path(base)/s).resolve()
def accession_key(s):return re.sub(r'\.\d+$','',s.strip().upper())
def fasta(p):
 header=None; parts=[]
 with Path(p).open(encoding='ascii') as f:
  for line in f:
   x=line.strip()
   if not x:continue
   if x.startswith('>'):
    if header is not None:raise Hold('multi-record FASTA; provide each complete native contig as a separate frozen row')
    header=x[1:]
   else:
    if header is None:raise Hold('sequence before FASTA header')
    parts.append(x.upper())
 if not header:raise Hold('missing FASTA header')
 seq=''.join(parts)
 if not seq:raise Hold('empty FASTA')
 return header,seq
def sequence_summary(p):
 _,s=fasta(p);rc=s.translate(str.maketrans('ACGT','TGCA'))[::-1]
 n=len(s);nw=0;largest_count=0;window_chars=0
 for L in POLICY['windows_bp']:
  if n<L:continue
  stride=int(L*.25); count=(n-L)//stride+1
  if (n-L)%stride:count+=1
  nw+=2*count
  largest_count=max(largest_count,count);window_chars+=2*count*L
 dense=largest_count*4096*4
 return {'length_bp':n,'sequence_sha256':digest(s),'reverse_complement_sha256':digest(rc),
   'fasta_sha256':sha(p),'non_acgt_count':sum(x not in 'ACGT' for x in s),
   'oriented_windows':nw,'raw_memmap_bytes':nw*5468*4,
   'full_legacy_raw_memmap_bytes':nw*(5468+3570)*4,
   'raw_memmap_execution_mode':'exact_two_arm',
   'largest_float32_scoring_block_bytes':dense,
   'planning_ram_bytes_heuristic_not_measured':8*dense+window_chars+nw*512,
   'gc_fraction':(s.count('G')+s.count('C'))/n}

def validate_manifest(candidates,truth,exclusions,analysis_kind="primary_art"):
 """Validate independently supplied evidence. Never call any detector here."""
 cbase=Path(candidates).resolve().parent;tbase=Path(truth).resolve().parent
 cc=read_tsv(candidates);tt=read_tsv(truth);xx=read_tsv(exclusions)
 hard_accessions={accession_key(x['value']) for x in xx if x['key_type'] in ('accession','native_parent_accession','contig_id')}
 hard_hashes={x['value'] for x in xx if x['key_type'] in ('sequence_sha256','reverse_complement_sha256')}
 known_families={x['value'] for x in xx if x['key_type']=='family_cluster'}
 tr={};errors=[];notes=[];obs=[];seen_ids=set();seen_seq=set();secondary=analysis_kind=='secondary_structured_rna'
 seen_targets=set()
 for t in tt:
  sid=t.get('sample_id','')
  target=(sid,t.get('target_id',sid))
  if target in seen_targets:errors.append({'sample_id':sid,'reason':'duplicate target_id within sample'})
  seen_targets.add(target)
  tr.setdefault(sid,[]).append(t)
 if not cc:errors.append({'sample_id':'','reason':'no externally eligible candidates supplied'})
 for c in cc:
  sid=c.get('sample_id',''); e=[];z={'sample_id':sid};sample_truth=tr.get(sid,[]);t=sample_truth[0] if sample_truth else {}
  if not re.fullmatch(r'[A-Za-z0-9_.-]+',sid):e.append('invalid sample_id')
  if sid in seen_ids:e.append('duplicate sample_id')
  seen_ids.add(sid)
  for k in [x for x in CANDIDATE_FIELDS if x!='accession_aliases']:
   if not c.get(k):e.append('missing candidate field '+k)
  for k in ['sample_id','label_class','label_scope','truth_basis','curator','truth_fixed_utc']:
   if not t.get(k):e.append('missing truth field '+k)
  if not secondary and c.get('exposure_status')!='UNEXPOSED':e.append('exposure_status must be independently reviewed UNEXPOSED')
  if c.get('eligibility_status') not in ({'SECONDARY_ONLY','CONFIRMED'} if secondary else {'CONFIRMED'}):e.append('eligibility_status incompatible with analysis kind')
  if secondary:notes.append({'sample_id':sid,'note':'secondary structure transfer only; exposure status '+c.get('exposure_status','UNKNOWN')+'; not confirmatory ART, not an ART negative'})
  if c.get('leakage_review_status')!='REVIEWED_SOURCE_AND_SEQUENCE_GROUPS':e.append('source/sequence/family leakage review unfinished')
  for key in ['accession','native_parent_accession']:
   if accession_key(c.get(key,'')) in hard_accessions:e.append('previously exposed '+key) if not secondary else notes.append({'sample_id':sid,'note':'previously exposed '+key+'; secondary/exposed status mandatory'})
  for alias in re.split(r'[;,|]',c.get('accession_aliases','')):
   if alias.strip() and accession_key(alias) in hard_accessions:
    if not secondary:e.append('previously exposed accession alias '+alias.strip())
    else:notes.append({'sample_id':sid,'note':'previously exposed alias '+alias.strip()+'; secondary only'})
  if c.get('family_cluster') in known_families:
   notes.append({'sample_id':sid,'note':'FAMILY_EXPOSED: family-disjoint generalization prohibited; independent accession/source evaluation only, subject to separate evidence review'})
  for pk,hk,base,row in [('leakage_review_path','leakage_review_sha256',cbase,c)]:
   try:
    p=resolved(row[pk],base)
    if sha(p)!=row[hk]:e.append(pk+' hash mismatch')
   except (KeyError,OSError) as ex:e.append(pk+' unavailable: '+str(ex))
  try:
   p=resolved(c['genome_path'],cbase);z.update(sequence_summary(p));z['genome_path']=str(p)
   if z['length_bp']<POLICY['minimum_sequence_length_bp']:e.append('native contig shorter than frozen minimum 4096 bp')
   if z['non_acgt_count']:e.append('ambiguous bases; no silent masking, deletion, or splitting permitted')
   for key in ['fasta_sha256','sequence_sha256']:
    if z[key]!=c.get(key):e.append(key+' mismatch')
   if z['sequence_sha256'] in hard_hashes or z['reverse_complement_sha256'] in hard_hashes:
    if not secondary:e.append('sequence or reverse complement previously exposed')
    else:notes.append({'sample_id':sid,'note':'exact native sequence previously exposed; secondary/exposed only'})
   if z['sequence_sha256'] in seen_seq or z['reverse_complement_sha256'] in seen_seq:e.append('duplicate sequence within new cohort')
   seen_seq.update([z['sequence_sha256'],z['reverse_complement_sha256']])
  except (KeyError,OSError,UnicodeError,Hold) as ex:e.append('FASTA: '+str(ex))
  cl=t.get('label_class');scope=t.get('label_scope')
  allowed={'other_structured_rna','background_unknown'} if secondary else {'positive','hard_neighbor','certified_negative','background_unknown'}
  if cl not in allowed:e.append('invalid label_class')
  basis={'other_structured_rna':'INDEPENDENT_AUTHOR_ANNOTATION','positive':'INDEPENDENT_AUTHOR_ANNOTATION','hard_neighbor':'INDEPENDENT_NON_ART_ANNOTATION',
         'certified_negative':'INDEPENDENT_NEGATIVE_CERTIFICATION','background_unknown':'NO_ART_NEGATIVE_CLAIM'}
  if cl in basis and t.get('truth_basis')!=basis[cl]:e.append('label lacks independent source definition; detector-derived truth forbidden')
  if cl!='background_unknown':
   for k in ['author_truth_url','truth_evidence_path','truth_evidence_sha256']:
    if not t.get(k):e.append('missing truth evidence '+k)
   try:
    if sha(resolved(t['truth_evidence_path'],tbase))!=t['truth_evidence_sha256']:e.append('truth evidence hash mismatch')
   except (KeyError,OSError) as ex:e.append('truth evidence unavailable: '+str(ex))
  try:
   d=dt.datetime.fromisoformat(t['truth_fixed_utc'].replace('Z','+00:00'))
   if not d.tzinfo or d>dt.datetime.now(dt.timezone.utc):e.append('truth_fixed_utc not a past timezone-aware timestamp')
  except (KeyError,ValueError):e.append('invalid truth_fixed_utc')
  if cl in ('positive','hard_neighbor','other_structured_rna'):
   if scope not in ({'author_defined_rna_interval'} if secondary else {'author_defined_system_interval'}):e.append('positive/hard_neighbor scope must be independent author-defined system interval')
   try:
    s=int(t['truth_start0']);en=int(t['truth_end0'])
    if not 0<=s<en<=z['length_bp']:e.append('invalid truth interval 0-based half-open')
   except (KeyError,ValueError):e.append('truth interval unresolved')
  if cl=='certified_negative' and scope!='whole_native_contig':e.append('certified-negative evidence must certify whole native contig')
  if cl=='background_unknown' and scope!='unlabelled_whole_native_contig':e.append('unknown background scope mismatch')
  for prefix in ['array','rt']:
   av=t.get(prefix+'_start0','');bv=t.get(prefix+'_end0','')
   if av or bv:
    try:
     a=int(av);b=int(bv)
     if not 0<=a<b<=z['length_bp']:e.append(prefix+' coordinates invalid')
    except (ValueError,KeyError):e.append(prefix+' coordinates incomplete')
  # Every extra independently annotated target is validated, not silently dropped.
  for extra in sample_truth[1:]:
   if extra.get('label_class')!=cl or extra.get('truth_basis')!=basis.get(cl):e.append('mixed target classes or non-independent truth within native sample')
   try:
    if not 0<=int(extra['truth_start0'])<int(extra['truth_end0'])<=z['length_bp']:e.append('invalid additional target interval')
    ep=resolved(extra['truth_evidence_path'],tbase)
    if sha(ep)!=extra['truth_evidence_sha256']:e.append('additional target evidence hash mismatch')
    ed=dt.datetime.fromisoformat(extra['truth_fixed_utc'].replace('Z','+00:00'))
    if not ed.tzinfo or ed>dt.datetime.now(dt.timezone.utc):e.append('invalid additional target timestamp')
   except (KeyError,ValueError,OSError):e.append('additional target evidence/coordinates incomplete')
  # Apply the complete scope and component-coordinate contract to every target.
  # Extra targets are never admitted with weaker evidence than the first row.
  for target in sample_truth:
   target_id=target.get('target_id',sid)
   if not target.get('target_id'):e.append('missing target_id for '+sid)
   for key in ['label_class','label_scope','truth_basis','curator','truth_fixed_utc']:
    if not target.get(key):e.append(target_id+' missing '+key)
   label=target.get('label_class')
   expected_scope={'positive':'author_defined_system_interval','hard_neighbor':'author_defined_system_interval',
     'other_structured_rna':'author_defined_rna_interval','certified_negative':'whole_native_contig',
     'background_unknown':'unlabelled_whole_native_contig'}.get(label)
   if target.get('label_scope')!=expected_scope:e.append(target_id+' label_scope mismatch')
   if label!='background_unknown':
    for key in ['author_truth_url','truth_evidence_path','truth_evidence_sha256']:
     if not target.get(key):e.append(target_id+' missing '+key)
   for prefix in ['array','rt']:
    av=target.get(prefix+'_start0','');bv=target.get(prefix+'_end0','')
    if label=='positive' and (not av or not bv):e.append(target_id+' ART-positive requires independent '+prefix+' coordinates')
    if av or bv:
     try:
      start=int(av);end=int(bv)
      if not 0<=start<end<=z['length_bp']:e.append(target_id+' invalid '+prefix+' coordinates')
      if label=='positive' and not int(target['truth_start0'])<=start<end<=int(target['truth_end0']):e.append(target_id+' '+prefix+' component outside author-defined system interval')
     except (ValueError,KeyError):e.append(target_id+' incomplete '+prefix+' coordinates')
  z['admitted']=not e;z['reasons']='; '.join(e);obs.append(z)
  errors.extend({'sample_id':sid,'reason':x} for x in e)
 extras=sorted(set(tr)-seen_ids)
 if extras:errors.append({'sample_id':','.join(extras),'reason':'truth rows absent from candidate cohort'})
 # Matched control membership is fixed before score, not selected after performance.
 groups={}
 for c in cc:groups.setdefault(c.get('match_group',''),set()).add((tr.get(c.get('sample_id',''),[{}])[0]).get('label_class',''))
 for group,classes in groups.items():
  if 'positive' in classes and not classes.intersection({'hard_neighbor','certified_negative'}):
   errors.append({'sample_id':group,'reason':'positive match_group lacks independently annotated matched hard-neighbor or certified-negative control'})
 if not secondary and not any(t.get('label_class')=='positive' for t in tt):errors.append({'sample_id':'','reason':'no independently annotated positive target'})
 # All samples sharing either source study, family or native parent must share one bootstrap unit.
 for key in ['native_parent_accession','family_cluster','source_study','match_group']:
  grouping={}
  for c in cc:grouping.setdefault(c.get(key,''),set()).add(c.get('leakage_unit',''))
  for value,units in grouping.items():
   if value and len(units)>1:errors.append({'sample_id':value,'reason':key+' split across leakage units; merge connected dependence components before freeze'})
 return cc,tt,obs,errors,notes

def admission(a):
 cc,tt,obs,errors,notes=validate_manifest(a.candidates,a.truth,a.exclusions,a.analysis_kind)
 out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
 write_tsv(out/'ADMISSION_INVENTORY.tsv',obs,['sample_id','admitted','length_bp','fasta_sha256','sequence_sha256','reverse_complement_sha256','non_acgt_count','gc_fraction','oriented_windows','raw_memmap_bytes','largest_float32_scoring_block_bytes','planning_ram_bytes_heuristic_not_measured','genome_path','reasons'])
 write_tsv(out/'ADMISSION_BLOCKERS.tsv',errors,['sample_id','reason']);write_tsv(out/'ADMISSION_NOTES.tsv',notes,['sample_id','note'])
 result={'status':'ELIGIBLE_FOR_PREFREEZE' if not errors else 'HOLD_EXTERNAL_ADMISSION',
  'analysis_kind':a.analysis_kind,'candidate_n':len(cc),'truth_n':len(tt),'blocker_n':len(errors),'scientific_result':False,
  'no_detector_used_for_labels':True,'family_disjoint_generalization_certified':False,
  'family_exposure_notes':[x for x in notes if 'FAMILY_EXPOSED' in x['note']], 'created_utc':now()}
 write_json(out/'ADMISSION_STATUS.json',result);print(json.dumps(result,indent=2));return 0 if not errors else 3

def verify_manifest(root,manifest):
 for r in read_tsv(manifest):
  p=Path(root)/r['path']
  if not p.is_file() or sha(p)!=r['sha256']:raise Hold('frozen file changed or missing: '+str(p))
def hashes(root,paths):return [{'path':str(p.relative_to(root)),'sha256':sha(p)} for p in sorted(paths) if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc']
def tool_record(cmd,flag):
 p=shutil.which(cmd) or (str(Path(cmd).resolve()) if Path(cmd).is_file() else None)
 if not p:raise Hold('required executable unavailable: '+cmd)
 q=subprocess.run([p,flag],capture_output=True,text=True,timeout=30)
 return {'path':str(Path(p).resolve()),'sha256':sha(p),'version_text':(q.stdout+q.stderr)[:1500]}
def find_rvt(project_root,explicit=None):
 if explicit:paths=[Path(explicit)]
 else:
  root=Path(project_root);patterns=[
   'ART_MATH_DISCOVERY_STAGE05D1_NOVELTY_FIRST_MINING*/references/myRT/RVT-All.hmm',
   'ART_MATH_DISCOVERY_STAGE05D4R5B_RT_FOCUSED_STRUCTURAL_PANEL*/references/RVT-All.hmm',
   'ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY*/references/RVT-All.hmm']
  paths=[p for pat in patterns for p in root.glob(pat)]
 for p in paths:
  if p.is_file() and sha(p)==RVT_SHA:return p.resolve()
 raise Hold('original RVT-All.hmm with SHA256 '+RVT_SHA+' not found; no replacement permitted')

def freeze(a):
 run=Path(a.run).resolve()
 if run.exists():raise Hold('run already exists; frozen runs are immutable, choose a new run directory')
 cc,tt,obs,errors,notes=validate_manifest(a.candidates,a.truth,a.exclusions,a.analysis_kind)
 if errors:raise Hold('external admission has '+str(len(errors))+' blockers; run admit for detailed report')
 if a.cpu<1:raise Hold('cpu must be positive')
 verify_manifest(HERE,HERE/'reference/VENDORED_SOURCE_SHA256.tsv')
 prod=tool_record(a.prodigal,'-v');hm=tool_record(a.hmmscan,'-h');press=tool_record(a.hmmpress,'-h')
 rvt=find_rvt(a.project_root,a.rvt)
 versions={x:importlib.metadata.version(x) for x in DEPENDENCIES}
 for name in DEPENDENCIES:importlib.import_module(name)
 run.mkdir(parents=True);f=run/'frozen';f.mkdir();sealed=run/'sealed';sealed.mkdir()
 shutil.copytree(HERE/'vendor',f/'vendor',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 shutil.copytree(HERE/'scripts',f/'driver',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 inp=f/'inputs';inp.mkdir();sc=[];mapping=[];ob={x['sample_id']:x for x in obs}
 for i,c in enumerate(sorted(cc,key=lambda x:x['sample_id']),1):
  oid='J2_'+format(i,'05d');z=ob[c['sample_id']];_,seq=fasta(z['genome_path']);dst=inp/(oid+'.frozen.fasta')
  dst.write_text('>'+oid+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
  sc.append({'opaque_id':oid,'fasta_path':str(dst.relative_to(run)),'fasta_sha256':sha(dst),
   'sequence_sha256':z['sequence_sha256'],'length_bp':z['length_bp'],'oriented_windows':z['oriented_windows'],
   'raw_memmap_bytes':z['full_legacy_raw_memmap_bytes'] if a.execution_mode=='full_legacy' else z['raw_memmap_bytes']})
  mapping.append({'opaque_id':oid,'sample_id':c['sample_id']})
 write_tsv(f/'SCORING_MANIFEST.tsv',sc);write_tsv(sealed/'IDENTITY_KEY.tsv',mapping)
 shutil.copyfile(a.truth,sealed/'TRUTH.tsv');shutil.copyfile(a.candidates,sealed/'CANDIDATES.tsv')
 shutil.copyfile(a.exclusions,sealed/'EXCLUSION_LEDGER.tsv')
 evidence_map=[]
 for typ,rows,base,pk,hk in [('truth',tt,Path(a.truth).resolve().parent,'truth_evidence_path','truth_evidence_sha256'),
   ('review',cc,Path(a.candidates).resolve().parent,'leakage_review_path','leakage_review_sha256')]:
  for row in rows:
   if row.get(pk):
    p=resolved(row[pk],base);h=sha(p);d=sealed/'evidence'/(typ+'_'+h+p.suffix);d.parent.mkdir(exist_ok=True)
    if not d.exists():shutil.copyfile(p,d)
    evidence_map.append({'kind':typ,'sample_id':row['sample_id'],'target_id':row.get('target_id',''),
      'original_path':str(p),'sealed_path':str(d.relative_to(run)),'sha256':h})
 write_tsv(sealed/'EVIDENCE_MAP.tsv',evidence_map)
 db=f/'db';db.mkdir();frozen_rvt=db/'RVT-All.hmm';shutil.copyfile(rvt,frozen_rvt)
 q=subprocess.run([press['path'],'-f',str(frozen_rvt)],capture_output=True,text=True)
 (db/'hmmpress.log').write_text(q.stdout+q.stderr)
 if q.returncode:raise Hold('hmmpress failed; incomplete run retained for diagnosis')
 for suffix in ['.h3f','.h3i','.h3m','.h3p']:
  if not Path(str(frozen_rvt)+suffix).is_file():raise Hold('hmmpress sidecar missing '+suffix)
 pol=dict(POLICY,analysis_kind=a.analysis_kind,cpu=a.cpu,per_case_timeout_seconds=a.timeout_seconds,python_version=platform.python_version(),package_versions=versions,
   execution_mode=a.execution_mode,exact_two_arm_adapter_sha256=sha(f/'driver/exact_two_arm_scan.py'))
 write_json(f/'POLICY.json',pol)
 write_json(f/'TOOLCHAIN.json',{'prodigal':prod,'hmmscan':hm,'hmmpress':press,
  'python':{'path':sys.executable,'sha256':sha(Path(sys.executable).resolve())},'rvt_source_sha256':sha(rvt)})
 write_tsv(f/'FROZEN_FILE_SHA256.tsv',hashes(run,[p for p in f.rglob('*') if p.name!='FROZEN_FILE_SHA256.tsv']))
 write_tsv(sealed/'SEALED_FILE_SHA256.tsv',hashes(run,[p for p in sealed.rglob('*') if p.name!='SEALED_FILE_SHA256.tsv']))
 gate={'status':'FROZEN_BEFORE_SCORING','frozen_utc':now(),'sample_n':len(sc),'cpu':a.cpu,
  'truth_seal_manifest_sha256':sha(sealed/'SEALED_FILE_SHA256.tsv'),'frozen_manifest_sha256':sha(f/'FROZEN_FILE_SHA256.tsv'),
  'scoring_reads_truth':False,'human_blinding_claimed':False,'fresh_external_performance_available':False,
  'admission_notes':notes,'analysis_kind':a.analysis_kind,'execution_mode':a.execution_mode,
  'exact_two_arm_adapter_sha256':pol['exact_two_arm_adapter_sha256'],'scientific_result':False}
 gate['family_disjoint_generalization_certified']=False
 write_json(run/'PREFREEZE.json',gate)
 print(json.dumps(gate,indent=2));return 0

def run_logged(cmd,log,timeout=0):
 log.parent.mkdir(parents=True,exist_ok=True)
 started=now();tic=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
 timer=Path('/usr/bin/time');timelog=log.with_suffix(log.suffix+'.time.txt')
 launched=([str(timer),'-v','-o',str(timelog),'--']+cmd) if timer.is_file() else cmd
 rc=None;timed_out=False
 child_usage=None
 with log.open('w') as f:
  proc=subprocess.Popen(launched,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  if hasattr(os,'wait4'):
   # wait4 provides actual per-command rusage even when GNU time is absent.
   # Popen.wait would consume this information, so reap this one child here.
   kill_deadline=None
   while True:
    pid,status,usage=os.wait4(proc.pid,os.WNOHANG)
    if pid:
     child_usage=usage;proc.returncode=os.waitstatus_to_exitcode(status);rc=proc.returncode;break
    if timeout and not timed_out and time.monotonic()-tic>=timeout:
     os.killpg(proc.pid,signal.SIGTERM);timed_out=True;kill_deadline=time.monotonic()+10
    if kill_deadline is not None and time.monotonic()>=kill_deadline:
     try:os.killpg(proc.pid,signal.SIGKILL)
     except ProcessLookupError:pass
     kill_deadline=None
    time.sleep(.1)
   if timed_out:
    try:os.killpg(proc.pid,signal.SIGKILL)
    except ProcessLookupError:pass
    rc=124
  else:
   try:rc=proc.wait(timeout=timeout or None)
   except subprocess.TimeoutExpired:
    os.killpg(proc.pid,signal.SIGTERM)
    try:proc.wait(timeout=10)
    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
    rc=124;timed_out=True
 after=resource.getrusage(resource.RUSAGE_CHILDREN)
 metrics={'started_utc':started,'finished_utc':now(),'command':cmd,'return_code':rc,'timed_out':timed_out,
   'wall_seconds':time.monotonic()-tic,'child_user_cpu_seconds':after.ru_utime-before.ru_utime,
   'child_system_cpu_seconds':after.ru_stime-before.ru_stime,
   'child_process_peak_rss_kib':None,'peak_rss_source':'not measured for this process if GNU time unavailable',
   'gnu_time_log':str(timelog) if timelog.exists() else None}
 if child_usage is not None:
  metrics['child_user_cpu_seconds']=child_usage.ru_utime
  metrics['child_system_cpu_seconds']=child_usage.ru_stime
  metrics['child_process_peak_rss_kib']=child_usage.ru_maxrss
  metrics['peak_rss_source']='Linux wait4 rusage; maximum resident set of the reaped child/descendants, not sum of concurrent worker RSS'
 if timelog.exists():
  match=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',timelog.read_text())
  if match:metrics['child_process_peak_rss_kib']=int(match.group(1));metrics['peak_rss_source']='GNU time -v; maximum child process RSS, not sum across simultaneous workers'
 write_json(log.with_suffix(log.suffix+'.resources.json'),metrics)
 return rc

def verify_frozen(run):
 gate=js(run/'PREFREEZE.json');f=run/'frozen'
 frozen_driver=f/'driver/j2_pipeline.py'
 if not frozen_driver.is_file() or sha(Path(__file__).resolve())!=sha(frozen_driver):raise Hold('current coordinator differs from the prefreeze copy; use the matching frozen driver')
 if sha(f/'FROZEN_FILE_SHA256.tsv')!=gate['frozen_manifest_sha256']:raise Hold('frozen manifest changed')
 verify_manifest(run,f/'FROZEN_FILE_SHA256.tsv')
 chain=js(f/'TOOLCHAIN.json')
 if sha(Path(sys.executable).resolve())!=chain['python']['sha256']:raise Hold('current Python executable differs from prefreeze interpreter')
 for x in chain.values():
  if isinstance(x,dict) and 'path' in x and sha(x['path'])!=x['sha256']:raise Hold('frozen executable changed: '+x['path'])
 pol=js(f/'POLICY.json')
 for name,version in pol['package_versions'].items():
  if importlib.metadata.version(name)!=version:raise Hold('frozen dependency version changed: '+name)
 return pol

def score(a):
 """Only frozen/ and scores/ are read; sealed/ is deliberately never opened."""
 run=Path(a.run).resolve();p=verify_frozen(run);out=run/'scores'
 if out.exists():raise Hold('scores directory already exists; do not overwrite or rerun after exposure')
 out.mkdir();f=run/'frozen';manifest=read_tsv(f/'SCORING_MANIFEST.tsv');statuses=[]
 adapter=f/'driver/score_adapter.py';cpu=p['cpu'];timeout=p['per_case_timeout_seconds']
 for row in manifest:
  cid=row['opaque_id'];cd=out/cid;cd.mkdir();arm=cd/'arm'
  shutil.copytree(f/'vendor/arm_m/config',arm/'config');shutil.copytree(f/'vendor/arm_m/frozen_sources',arm/'frozen_sources')
  (arm/'inputs').mkdir();shutil.copyfile(run/row['fasta_path'],arm/'inputs'/(cid+'.frozen.fasta'))
  if p['execution_mode']=='full_legacy':
   cmd=[sys.executable,str(f/'vendor/arm_m/scripts/02_primary_scan.py')]
  elif p['execution_mode']=='exact_two_arm':
   cmd=[sys.executable,str(f/'driver/exact_two_arm_scan.py'),'--vendor',str(f/'vendor')]
  else:raise Hold('unknown frozen execution mode')
  rc=run_logged(cmd+['--root',str(arm),'--accession',cid,'--threads',str(cpu),'--keep-raw','0'],cd/'arm_scan.log',timeout)
  if rc==0:
   rc=run_logged([sys.executable,str(adapter),'--run',str(run),'--id',cid,'--phase','board'],cd/'board.log',timeout)
  arm_ok=rc==0
  statuses.append({'opaque_id':cid,'method':METHODS[0],'status':'SUCCESS' if arm_ok else 'FAILED',
   'exit_code':rc,'output_path':str((cd/(METHODS[0]+'.tsv')).relative_to(run)),'failure_reason':'' if arm_ok else 'raw-DNA scan or frozen board build failed; see logs'})
  for route,prefix in [('prodigal','PRODIGAL_RVT'),('sixframe','SIXFRAME_RVT')]:
   brc=run_logged([sys.executable,str(adapter),'--run',str(run),'--id',cid,'--phase',route],cd/(route+'.log'),timeout)
   route_status=cd/'independent_rt'/route/'ROUTE_STATUS.json'
   rec_available=js(route_status).get('recurrence_available',False) if brc==0 and route_status.is_file() else False
   for suffix in ['_PLUS_FROZEN_REC','_ONLY']:
    method=prefix+suffix;rec_missing=suffix=='_PLUS_FROZEN_REC' and not rec_available;ok=brc==0 and not rec_missing
    statuses.append({'opaque_id':cid,'method':method,'status':'SUCCESS' if ok else 'FAILED',
     'exit_code':125 if rec_missing and brc==0 else brc,'output_path':str((cd/(method+'.tsv')).relative_to(run)),
     'failure_reason':'' if ok else ('DEPENDENCY_FAILURE: frozen recurrence unavailable; RT-only is still independently attempted' if rec_missing and brc==0 else route+' RT engine failed; see independent route log')})
  write_tsv(out/'METHOD_STATUS.tsv',statuses)
  print(json.dumps({'opaque_id':cid,'methods':{x['method']:x['status'] for x in statuses if x['opaque_id']==cid}}),flush=True)
 # Recheck mutable environment and all frozen scoring resources after execution.
 verify_frozen(run)
 write_tsv(out/'SCORE_FILE_SHA256.tsv',hashes(run,[x for x in out.rglob('*') if x.name!='SCORE_FILE_SHA256.tsv']))
 gate={'status':'SCORES_LOCKED_ALL_METHOD_ATTEMPTS_RECORDED','score_frozen_utc':now(),
  'sample_n':len(manifest),'method_n':len(METHODS),'expected_attempt_n':len(manifest)*len(METHODS),
  'recorded_attempt_n':len(statuses),'success_n':sum(x['status']=='SUCCESS' for x in statuses),
  'failed_n':sum(x['status']!='SUCCESS' for x in statuses),'score_manifest_sha256':sha(out/'SCORE_FILE_SHA256.tsv'),
  'truth_read_for_scoring':False,'scientific_performance_not_yet_evaluated':True}
 write_json(run/'SCORE_FREEZE.json',gate);print(json.dumps(gate,indent=2));return 0

def overlap(a,b,c,d):return a<d and c<b
def union_length(xs):
 total=0;last=None
 for a,b in sorted(xs):
  if last is None:last=[a,b]
  elif a<=last[1]:last[1]=max(last[1],b)
  else:total+=last[1]-last[0];last=[a,b]
 return total+(last[1]-last[0] if last else 0)
def endpoint(rows,start,end):
 ranks=[int(r['rank']) for r in rows if overlap(int(r['start0']),int(r['end0']),start,end)]
 rank=min(ranks) if ranks else 0
 return {'rank':rank or '', 'recovered':int(rank>0),'hit1':int(rank==1),
  'hit3':int(0<rank<=3),'hit5':int(0<rank<=5),'mrr':1/rank if rank else 0.0}
def quantile(xs,q):
 x=sorted(xs);z=(len(x)-1)*q;i=int(z);return x[i]+(x[min(i+1,len(x)-1)]-x[i])*(z-i)

def summarize(rows,pol):
 pos=[r for r in rows if r['label_class'] in ('positive','other_structured_rna')];result=[]
 for method in METHODS:
  mm=[r for r in pos if r['method']==method];n=len(mm)
  result.append({'method':method,'role':'primary' if method in PRIMARY else 'secondary',
   'target_denominator':n,'successful_target_n':sum(r['status']=='SUCCESS' for r in mm),
   **{x:sum(float(r[x]) for r in mm)/n if n else None for x in ['recovered','hit1','hit3','hit5','mrr']}})
 units=sorted(set(r['leakage_unit'] for r in pos));ci=[]
 if len(units)>=pol['bootstrap_minimum_positive_units']:
  # Same sampled dependence units in all methods and paired contrasts.
  by={u:[r for r in pos if r['leakage_unit']==u] for u in units};rng=random.Random(pol['bootstrap_seed']);vals={}
  for _ in range(pol['bootstrap_replicates']):
   sampled=[r for u in rng.choices(units,k=len(units)) for r in by[u]];means={}
   for m in METHODS:
    z=[r for r in sampled if r['method']==m]
    for met in ['recovered','hit1','hit3','hit5','mrr']:
     means[(m,met)]=sum(float(r[met]) for r in z)/len(z)
     vals.setdefault((m,met),[]).append(means[(m,met)])
   for m in PRIMARY[1:]:
    for met in ['recovered','hit1','hit3','hit5','mrr']:
     vals.setdefault(('ARM_M_FULL_FROZEN minus '+m,met),[]).append(means[(PRIMARY[0],met)]-means[(m,met)])
  for (m,met),v in vals.items():ci.append({'method_or_contrast':m,'metric':met,'lower_0.025':quantile(v,.025),'upper_0.975':quantile(v,.975),'independent_positive_units':len(units),'replicates':pol['bootstrap_replicates'],'seed':pol['bootstrap_seed']})
 return result,ci,len(units)

def evaluate(a):
 run=Path(a.run).resolve();p=verify_frozen(run);out=run/'evaluation'
 if out.exists():raise Hold('evaluation already exists; keep locked result immutable')
 fg=js(run/'SCORE_FREEZE.json');scores=run/'scores'
 if sha(scores/'SCORE_FILE_SHA256.tsv')!=fg['score_manifest_sha256']:raise Hold('score freeze manifest changed')
 verify_manifest(run,scores/'SCORE_FILE_SHA256.tsv')
 ss=read_tsv(scores/'METHOD_STATUS.tsv');sm=read_tsv(run/'frozen/SCORING_MANIFEST.tsv')
 expected={(r['opaque_id'],m) for r in sm for m in METHODS};observed=[(r['opaque_id'],r['method']) for r in ss]
 if len(observed)!=len(expected) or set(observed)!=expected:raise Hold('incomplete method attempts; truth remains sealed')
 for r in ss:
  if r['status'] not in ('SUCCESS','FAILED'):raise Hold('unresolved method status; truth remains sealed')
  if r['status']=='SUCCESS' and not (run/r['output_path']).is_file():raise Hold('successful score output missing; truth remains sealed')
 # Only now is the sealed side opened. Hash verification detects altered labels.
 seal=run/'sealed';pref=js(run/'PREFREEZE.json')
 if sha(seal/'SEALED_FILE_SHA256.tsv')!=pref['truth_seal_manifest_sha256']:raise Hold('sealed truth manifest changed')
 verify_manifest(run,seal/'SEALED_FILE_SHA256.tsv')
 key={r['opaque_id']:r['sample_id'] for r in read_tsv(seal/'IDENTITY_KEY.tsv')}
 truth={}
 for r in read_tsv(seal/'TRUTH.tsv'):truth.setdefault(r['sample_id'],[]).append(r)
 candidates={r['sample_id']:r for r in read_tsv(seal/'CANDIDATES.tsv')}
 lengths={r['opaque_id']:int(r['length_bp']) for r in sm};evaluated=[];burden=[];components=[]
 for status in ss:
  oid=status['opaque_id'];sid=key[oid];t=truth[sid][0];c=candidates[sid];m=status['method'];n=lengths[oid]
  arr=read_tsv(run/status['output_path']) if status['status']=='SUCCESS' else []
  ranks=[int(r['rank']) for r in arr]
  if sorted(ranks)!=list(range(1,len(arr)+1)):raise Hold('malformed candidate ranks '+oid+'/'+m)
  for r in arr:
   if not 0<=int(r['start0'])<int(r['end0'])<=n:raise Hold('invalid candidate interval '+oid+'/'+m)
  base={'sample_id':sid,'opaque_id':oid,'method':m,'status':status['status'],'label_class':t['label_class'],
    'match_group':c['match_group'],'leakage_unit':c['leakage_unit'],'source_study':c['source_study'],'family_cluster':c['family_cluster']}
  measured=status['status']=='SUCCESS';iv=[(int(r['start0']),int(r['end0'])) for r in arr]
  er=[(max(0,a-10000),min(n,b+10000)) for a,b in iv] if m!=METHODS[0] else iv
  burden.append(dict(base,length_bp=n,candidate_n=len(arr) if measured else '',
   candidates_per_mb=len(arr)/(n/1e6) if measured else '',nominated_span_fraction=union_length(iv)/n if measured else '',
   rt_neighborhood_exposure_fraction=union_length(er)/n if measured and m!=METHODS[0] else '',
   interpretation='candidate burden; not a population false-positive rate or precision'))
  for t in truth[sid]:
   if t['label_class'] in ('positive','other_structured_rna'):
    tb=dict(base,target_id=t.get('target_id',sid))
    z=endpoint(arr,int(t['truth_start0']),int(t['truth_end0']));evaluated.append(dict(tb,**z))
    for target in ['array','rt']:
     if t.get(target+'_start0') and t.get(target+'_end0'):
      comp=endpoint(arr,int(t[target+'_start0']),int(t[target+'_end0']))
      components.append(dict(tb,target=target,**comp))
    if m!=METHODS[0] and t.get('array_start0') and t.get('array_end0'):
     expanded=[dict(r,start0=max(0,int(r['start0'])-10000),end0=min(n,int(r['end0'])+10000)) for r in arr]
     components.append(dict(tb,target='array_within_rt_10kb_neighborhood',**endpoint(expanded,int(t['array_start0']),int(t['array_end0']))))
 summaries,ci,units=summarize(evaluated,p);out.mkdir()
 per_genome=[]
 for sid in sorted(candidates):
  for m in METHODS:
   gg=[r for r in evaluated if r['sample_id']==sid and r['method']==m]
   if not gg:continue
   per_genome.append({'sample_id':sid,'method':m,'leakage_unit':gg[0]['leakage_unit'],
     'target_n':len(gg),'any_target_recovered':int(any(r['recovered'] for r in gg)),
     **{'mean_target_'+k:sum(float(r[k]) for r in gg)/len(gg) for k in ['recovered','hit1','hit3','hit5','mrr']}})
 genome_summary=[]
 for m in METHODS:
  gg=[r for r in per_genome if r['method']==m]
  genome_summary.append({'method':m,'native_genome_denominator':len(gg),
   'genomes_with_any_target_recovered':sum(r['any_target_recovered'] for r in gg),
   **{'genome_mean_'+k:sum(float(r['mean_target_'+k]) for r in gg)/len(gg) if gg else '' for k in ['recovered','hit1','hit3','hit5','mrr']}})
 write_tsv(out/'TARGET_METRICS.tsv',evaluated)
 write_tsv(out/'ALL_SAMPLE_CANDIDATE_BURDEN.tsv',burden)
 write_tsv(out/'COMPONENT_LOCALIZATION.tsv',components)
 write_tsv(out/'TARGET_SUMMARY.tsv',summaries)
 write_tsv(out/'PER_GENOME_TARGET_METRICS.tsv',per_genome)
 write_tsv(out/'GENOME_AGGREGATE_SUMMARY.tsv',genome_summary)
 write_tsv(out/'CLUSTER_BOOTSTRAP.tsv',ci,['method_or_contrast','metric','lower_0.025','upper_0.975','independent_positive_units','replicates','seed'])
 status_label='POSTFREEZE_SECONDARY_STRUCTURED_RNA_EVALUATION_COMPLETE' if p['analysis_kind']=='secondary_structured_rna' else 'POSTFREEZE_EXTERNAL_EVALUATION_COMPLETE'
 report={'status':status_label+('_WITH_RECORDED_FAILURES' if fg['failed_n'] else ''),
  'evaluated_utc':now(),'analysis_kind':p['analysis_kind'],'target_n':len(evaluated)//len(METHODS),'sample_n':len(sm),'method_attempt_n':len(ss),'failed_attempt_n':fg['failed_n'],
  'independent_positive_leakage_units':units,'bootstrap_status':'REPORTED' if ci else 'DESCRIPTIVE_ONLY_TOO_FEW_INDEPENDENT_POSITIVE_UNITS',
  'truth_read_only_after_all_score_attempts_locked':True,'population_specificity_claimed':False,
  'functional_validation_claimed':False,'score_manifest_sha256':fg['score_manifest_sha256']}
 if p['analysis_kind']=='secondary_structured_rna':
  report['interpretation']='General structured-RNA localization challenge. These targets need not contain RT; absent RT-route hits cannot be interpreted as poor ART detection. No ART superiority claim from this secondary cohort.'
  report['art_detection_superiority_claim_permitted']=False
  report['formal_art_method_comparison_requires_separate_primary_art_cohort']=True
 write_json(out/'EVALUATION_STATUS.json',report);print(json.dumps(report,indent=2));return 0

def preflight(a):
 out=Path(a.output);out.mkdir(parents=True,exist_ok=True);problems=[];tools={}
 try:verify_manifest(HERE,HERE/'reference/VENDORED_SOURCE_SHA256.tsv')
 except Hold as ex:problems.append(str(ex))
 for x,arg in [('prodigal','-v'),('hmmscan','-h'),('hmmpress','-h')]:
  try:tools[x]=tool_record(getattr(a,x),arg)
  except (Hold,subprocess.SubprocessError) as ex:problems.append(str(ex))
 versions={}
 for x in DEPENDENCIES:
  try:
   versions[x]=importlib.metadata.version(x);importlib.import_module(x)
  except (importlib.metadata.PackageNotFoundError,ImportError) as ex:problems.append('Python dependency missing or not importable: '+x+'; '+str(ex))
 try:rv=find_rvt(a.project_root,a.rvt);tools['rvt']={'path':str(rv),'sha256':sha(rv)}
 except Hold as ex:problems.append(str(ex))
 memory={}
 for filename in ['/sys/fs/cgroup/memory.max','/sys/fs/cgroup/memory.current']:
  try:memory[Path(filename).name]=Path(filename).read_text().strip()
  except OSError:pass
 report={'status':'READY_FOR_ADMISSION_AND_FREEZE' if not problems else 'HOLD_RUNTIME_PREREQUISITES',
  'tools':tools,'python':platform.python_version(),'package_versions':versions,'blockers':problems,
  'memory_cgroup':memory,'output_filesystem_free_bytes':shutil.disk_usage(out).free,
  'resource_note':'admission estimates dense blocks and planning RAM; actual wait4/GNU time peak child RSS is recorded per stage; worker count is a memory as well as CPU choice',
  'external_scientific_truth_not_checked_here':True,'created_utc':now()}
 write_json(out/'PREFLIGHT.json',report);print(json.dumps(report,indent=2));return 0 if not problems else 3

def main():
 p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='command',required=True)
 for cmd in ['admit','freeze']:
  q=sp.add_parser(cmd);q.add_argument('--candidates',required=True);q.add_argument('--truth',required=True)
  q.add_argument('--exclusions',default=str(HERE/'reference/EXCLUSION_LEDGER.tsv'));q.add_argument('--analysis-kind',choices=['primary_art','secondary_structured_rna'],default='primary_art')
  if cmd=='admit':q.add_argument('--output',required=True)
  else:
   q.add_argument('--run',required=True);q.add_argument('--cpu',type=int,default=16);q.add_argument('--timeout-seconds',type=int,default=0)
   q.add_argument('--execution-mode',choices=['exact_two_arm','full_legacy'],default='exact_two_arm')
 for cmd in ['freeze','preflight']:
  q=sp.choices.get(cmd) or sp.add_parser(cmd)
  q.add_argument('--project-root',required=True);q.add_argument('--rvt')
  q.add_argument('--prodigal',default='prodigal');q.add_argument('--hmmscan',default='hmmscan');q.add_argument('--hmmpress',default='hmmpress')
  if cmd=='preflight':q.add_argument('--output',required=True)
 for cmd in ['score','evaluate']:
  q=sp.add_parser(cmd);q.add_argument('--run',required=True)
 a=p.parse_args()
 try:return {'admit':admission,'freeze':freeze,'score':score,'evaluate':evaluate,'preflight':preflight}[a.command](a)
 except (Hold,OSError,ValueError,KeyError,ImportError,importlib.metadata.PackageNotFoundError) as ex:
  print(json.dumps({'status':'HOLD','stage':a.command,'reason':str(ex),'scientific_success_claimed':False},indent=2),file=sys.stderr);return 3
if __name__=='__main__':sys.exit(main())
