#!/usr/bin/env python3
"""Real RNAfold only. Null fractions are descriptive, not calibrated p-values."""
from pathlib import Path
import argparse,collections,concurrent.futures,json,math,os,random,re,shutil,statistics,subprocess,sys
from common import *

def dinucs(s):return collections.Counter(zip(s,s[1:]))
def shuffled_dinucs(s,seed):
 if len(s)<2:return s
 rng=random.Random(seed);adj=collections.defaultdict(list)
 for x,y in zip(s,s[1:]):adj[x].append(y)
 for k in sorted(adj):rng.shuffle(adj[k])
 stack=[s[0]];trail=[]
 while stack:
  k=stack[-1]
  if adj[k]:stack.append(adj[k].pop())
  else:trail.append(stack.pop())
 out=''.join(reversed(trail))
 if len(out)!=len(s) or out[0]!=s[0] or out[-1]!=s[-1] or dinucs(out)!=dinucs(s):raise ValueError('dinucleotide shuffle failed')
 return out

def parse_output(text,expected):
 rows={};key=None;seq=None
 for ln in text.splitlines():
  ln=ln.strip()
  if not ln:continue
  if ln.startswith('>'):
   key=ln[1:].split()[0];seq=None
   if key not in expected or key in rows:raise ValueError('unknown/duplicate RNAfold ID '+key)
  elif key and seq is None:
   seq=ln.upper().replace('T','U')
   if seq!=expected[key]:raise ValueError('RNAfold sequence mismatch '+key)
  elif key:
   m=re.fullmatch(r'([.()]+)\s+\(\s*([-+]?\d+(?:\.\d+)?)\s*\)',ln)
   if not m:raise ValueError('unexpected RNAfold structure '+ln[:120])
   db,en=m.group(1),float(m.group(2));level=0
   if len(db)!=len(seq) or not math.isfinite(en):raise ValueError('invalid fold')
   for c in db:
    level+=(c=='(')-(c==')')
    if level<0:raise ValueError('unbalanced fold')
   if level:raise ValueError('unbalanced fold')
   rows[key]={'structure':db,'MFE_kcal_mol':en};key=None;seq=None
  else:raise ValueError('unparsed RNAfold line '+ln[:120])
 if set(rows)!=set(expected):raise ValueError('RNAfold incomplete results')
 return rows

def fold_all(root,jobs):
 root=Path(root);inp=root/'precomputed';out=root/'outputs';work=root/'work/fold';work.mkdir(parents=True,exist_ok=True);out.mkdir(exist_ok=True)
 policy=read_json(root/'config/POLICY.json');cfg=policy['RNAfold'];fr=read_json(inp/'INPUT_FREEZE.json')
 if sha_file(inp/'RNA_TARGETS.fasta')!=fr['RNA_targets_sha256'] or sha_file(inp/'RNA_OCCURRENCE_REGISTRY.tsv')!=fr['RNA_occurrences_sha256'] or sha_file(root/'config/POLICY.json')!=fr['policy_sha256']:raise ValueError('frozen inputs changed')
 exe=shutil.which('RNAfold')
 if not exe:raise RuntimeError('RNAfold is missing; activate art_rna_review (not the BLAST-only environment)')
 ver=subprocess.run([exe,'--version'],text=True,capture_output=True,check=True,timeout=20).stdout.strip()
 if not re.search(r'\b'+re.escape(cfg['required_version'])+r'\b',ver):raise RuntimeError('require RNAfold '+cfg['required_version']+'; observed '+ver)
 cmd=[str(Path(exe).resolve())]+cfg['arguments']
 run={'command':cmd,'version':ver,'binary_sha256':sha_file(exe),'policy_sha256':sha_file(root/'config/POLICY.json'),'targets_sha256':fr['RNA_targets_sha256'],'script_sha256':sha_file(__file__),'python_version':sys.version,'descriptive_not_calibrated_p':True}
 lock=out/'RNA_RUN_LOCK.json'
 if lock.exists() and read_json(lock)!=run:raise RuntimeError('run signature changed; do not mix runs')
 write_json(lock,run);targets=fasta(inp/'RNA_TARGETS.fasta');N=int(cfg['null_replicates'])
 def one(k):
  d=work/k;d.mkdir(exist_ok=True);cache=d/'result.json';seq=targets[k]
  seedprefix=cfg['seed_namespace']+'|'+k
  seqs={k+'__OBS':seq}
  for i in range(N):seqs[k+'__N%03d'%(i+1)]=shuffled_dinucs(seq,int(sha_bytes((seedprefix+'|'+str(i)).encode()),16))
  valid=False
  if cache.exists():
   r=read_json(cache);meta=r['cache_integrity']
   valid=meta['run_lock_sha256']==sha_file(lock) and meta['RNA_sha256']==sha_text(seq) and (d/'stdout.txt').is_file() and sha_file(d/'stdout.txt')==meta['stdout_sha256'] and (d/'input.fasta').is_file() and sha_file(d/'input.fasta')==meta['input_sha256']
   if not valid:raise RuntimeError('damaged completed target cache '+k)
   folds=parse_output((d/'stdout.txt').read_text(),seqs)
  else:
   write_fasta(d/'input.fasta',seqs)
   cp=subprocess.run(cmd,input=(d/'input.fasta').read_text(),capture_output=True,text=True,timeout=int(cfg['target_timeout_seconds']),cwd=d)
   (d/'stderr.txt').write_text(cp.stderr)
   if cp.returncode:raise RuntimeError('RNAfold failed '+k+': '+cp.stderr[:300])
   (d/'stdout.txt').write_text(cp.stdout);folds=parse_output(cp.stdout,seqs)
  obs=folds[k+'__OBS'];null=[folds[k+'__N%03d'%(i+1)]['MFE_kcal_mol'] for i in range(N)];oe=obs['MFE_kcal_mol'];mean=statistics.mean(null);sd=statistics.pstdev(null)
  r={'target_id':k,'RNA_length_nt':len(seq),'RNA_sha256':sha_text(seq),'MFE_kcal_mol':oe,'MFE_per_nt':oe/len(seq),'dot_bracket':obs['structure'],'base_pairs':obs['structure'].count('('),'null_n':N,'null_distinct_sequence_n':len(set(list(seqs.values())[1:])),'null_equal_to_observed_n':sum(x==seq for x in list(seqs.values())[1:]),'null_mean_MFE':mean,'null_sd_MFE':sd,'descriptive_standardized_energy':(oe-mean)/sd if sd>0 else None,'descriptive_lower_tail_fraction':(1+sum(x<=oe for x in null))/(1+N),'calibrated_p_value':False,'FDR_or_positive_call':False,'cache_integrity':{'run_lock_sha256':sha_file(lock),'RNA_sha256':sha_text(seq),'stdout_sha256':sha_file(d/'stdout.txt'),'input_sha256':sha_file(d/'input.fasta')}}
  if not valid:write_json(cache,r)
  return r
 results=[];errors=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
  fut={ex.submit(one,k):k for k in sorted(targets)}
  for f in concurrent.futures.as_completed(fut):
   try:
    r=f.result();results.append(r)
    if len(results)%25==0 or len(results)==len(targets):print('[RNAfold] %d/%d completed'%(len(results),len(targets)),flush=True)
   except Exception as e:errors.append({'target_id':fut[f],'error':type(e).__name__+':'+str(e)});print('[HOLD]',errors[-1],flush=True)
 results.sort(key=lambda r:r['target_id'])
 fields=[k for k in results[0] if k!='cache_integrity'] if results else ['target_id']
 write_tsv(out/'RNA_TARGET_FOLD_RESULTS.tsv',[{k:r[k] for k in fields} for r in results],fields)
 maps={r['target_id']:r for r in results};occ=read_tsv(inp/'RNA_OCCURRENCE_REGISTRY.tsv');oc=[]
 for r in occ:
  x=maps.get(r['target_id']);oc.append({**r,'fold_status':'COMPUTED' if x else 'MISSING','MFE_kcal_mol':x['MFE_kcal_mol'] if x else '', 'descriptive_lower_tail_fraction':x['descriptive_lower_tail_fraction'] if x else '', 'ART_positive_call':False})
 write_tsv(out/'ALL_RNA_OCCURRENCE_RESULTS.tsv',oc)
 write_json(out/'RNA_FOLD_GATE.json',{'gate':'PASS_ALL_FROZEN_RNA_MFE_RESULTS_WITH_SCOPE_LIMITS' if not errors and len(results)==len(targets) else 'HOLD_RNA_EXECUTION_INCOMPLETE','gate_pass':not errors and len(results)==len(targets),'n_expected_targets':len(targets),'n_completed_targets':len(results),'n_occurrences':len(occ),'null_replicates':N,'actual_RNAfold_version':ver,'RNAfold_executed':True,'errors':errors,'full_author_catalogue_comparison_complete':False,'RNA_expression_or_function_proven':False,'novel_ART_confirmed':0,'wet_lab_ready_declared':0,'UG27_touched':False})
 if errors:raise RuntimeError('RNA target errors: '+str(len(errors)))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--root',required=True);a.add_argument('--jobs',type=int,default=4);x=a.parse_args()
 if x.jobs<1:raise SystemExit('--jobs must be >= 1')
 try:fold_all(x.root,x.jobs)
 except Exception as e:
  write_json(Path(x.root)/'outputs/RNA_LAST_ERROR.json',{'gate':'HOLD_RNA_RUNTIME','error':type(e).__name__+': '+str(e),'no_biological_negative_call':True});raise
