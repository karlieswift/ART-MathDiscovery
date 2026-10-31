#!/usr/bin/env python3
"""Real full-genome CPU regression on previously exposed SA1 (260727 bp).

Recomputes the complete formal ACNV and recurrence features, then runs actual
Prodigal/HMMER and checks genomic boards against supplied historical evidence.
Optional full_legacy also computes the unused IRR branch. Not a new external
benchmark, not independent discovery, not a replacement for new-data admission.
"""
from __future__ import annotations
import argparse,json,math,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from j2_pipeline import *

def prepare(out,cpu,execution_mode):
 if out.exists():raise Hold('historical output already exists; choose a new output or use --phase finish for an already completed raw scan')
 ref=HERE/'reference/historical_SA1'
 for r in read_tsv(ref/'SOURCE_MANIFEST.tsv'):
  if sha(ref/r['filename'])!=r['sha256']:raise Hold('historical fixture changed '+r['filename'])
 verify_manifest(HERE,HERE/'reference/VENDORED_SOURCE_SHA256.tsv')
 out.mkdir(parents=True);shutil.copytree(HERE/'vendor',out/'frozen/vendor',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 fa=ref/'SA1__MW218148.1.frozen.fasta';header,seq=fasta(fa);cid='HIST_SA1';cd=out/'scores'/cid;arm=cd/'arm'
 shutil.copytree(HERE/'vendor/arm_m/config',arm/'config');shutil.copytree(HERE/'vendor/arm_m/frozen_sources',arm/'frozen_sources')
 (arm/'inputs').mkdir();inputfa=arm/'inputs'/(cid+'.frozen.fasta')
 inputfa.write_text('>'+cid+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
 write_tsv(out/'frozen/SCORING_MANIFEST.tsv',[{'opaque_id':cid,'fasta_path':str(inputfa.relative_to(out)),'length_bp':len(seq)}])
 write_json(out/'frozen/POLICY.json',dict(POLICY,cpu=cpu,execution_mode=execution_mode,
   exact_two_arm_adapter_sha256=sha(HERE/'scripts/exact_two_arm_scan.py'),analysis_kind='historical_regression_not_new_external',
   python_version=platform.python_version(),package_versions={x:importlib.metadata.version(x) for x in DEPENDENCIES}))
 write_json(out/'HISTORICAL_PROVENANCE.json',{'status':'RUNNING_TRUE_HISTORICAL_FULL_GENOME_REGRESSION',
  'original_fasta':str(fa),'original_fasta_sha256':sha(fa),'sequence_sha256':digest(seq),'original_header':header,
  'length_bp':len(seq),'new_external_validation':False,'truth_known_historically':True,
  'execution_mode':execution_mode,
  'raw_bytes_estimate':sequence_summary(fa)['raw_memmap_bytes' if execution_mode=='exact_two_arm' else 'full_legacy_raw_memmap_bytes']})

def finish(a,out):
 ref=HERE/'reference/historical_SA1';cid='HIST_SA1';cd=out/'scores'/cid;f=out/'frozen'
 for r in read_tsv(ref/'SOURCE_MANIFEST.tsv'):
  if sha(ref/r['filename'])!=r['sha256']:raise Hold('historical reference changed '+r['filename'])
 provenance=js(out/'HISTORICAL_PROVENANCE.json')
 actual_fa=cd/'arm/inputs'/f'{cid}.frozen.fasta'
 _,seq=fasta(actual_fa);_,expected_seq=fasta(ref/'SA1__MW218148.1.frozen.fasta')
 if seq!=expected_seq:raise Hold('historical input no longer equals complete SA1 native genome')
 q=cd/'arm/outputs'/cid/'STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz'
 if not q.is_file():raise Hold('raw historical scan incomplete')
 adapter=HERE/'scripts/score_adapter.py'
 rc=run_logged([sys.executable,str(adapter),'--run',str(out),'--id',cid,'--phase','board'],out/'board.log')
 if rc:raise Hold('historical frozen board adapter failed')
 if not a.rvt:raise Hold('historical baseline requires --rvt original hash-matched HMM')
 if sha(a.rvt)!=RVT_SHA:raise Hold('historical HMM fingerprint mismatch')
 tool={'prodigal':tool_record(a.prodigal,'-v'),'hmmscan':tool_record(a.hmmscan,'-h'),'hmmpress':tool_record(a.hmmpress,'-h')}
 write_json(f/'TOOLCHAIN.json',tool);(f/'db').mkdir(exist_ok=True);shutil.copyfile(a.rvt,f/'db/RVT-All.hmm')
 rc=run_logged([tool['hmmpress']['path'],'-f',str(f/'db/RVT-All.hmm')],out/'hmmpress.log')
 if rc:raise Hold('hmmpress failed')
 rc=run_logged([sys.executable,str(adapter),'--run',str(out),'--id',cid,'--phase','baseline'],out/'baseline.log')
 if rc:raise Hold('historical frozen baseline adapter failed')
 # Compare the new independent-route orchestration against the original paired
 # score_case wrapper, using actual tool executions and identical full DNA.
 paired_outputs={m:read_tsv(cd/(m+'.tsv')) for m in METHODS[1:]}
 route_checks=[]
 for route,prefix in [('prodigal','PRODIGAL_RVT'),('sixframe','SIXFRAME_RVT')]:
  rc=run_logged([sys.executable,str(adapter),'--run',str(out),'--id',cid,'--phase',route],out/(route+'_independent.log'))
  if rc:raise Hold('independent historical '+route+' adapter failed')
  for suffix in ['_PLUS_FROZEN_REC','_ONLY']:
   m=prefix+suffix;route_checks.append({'method':m,'condition':'recurrence_available','match_original_score_case':read_tsv(cd/(m+'.tsv'))==paired_outputs[m]})
 # A genuinely unavailable recurrence input must not suppress the two RT-only
 # computations. This is a dependency test on historical DNA, never a new case.
 dep=out/'missing_recurrence_dependency_test';df=dep/'frozen';df.mkdir(parents=True,exist_ok=True)
 shutil.copytree(f/'vendor',df/'vendor',dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 shutil.copytree(f/'db',df/'db',dirs_exist_ok=True)
 shutil.copyfile(f/'POLICY.json',df/'POLICY.json');shutil.copyfile(f/'TOOLCHAIN.json',df/'TOOLCHAIN.json')
 write_tsv(df/'SCORING_MANIFEST.tsv',[{'opaque_id':cid,'fasta_path':str(actual_fa),'length_bp':len(seq)}])
 for route,prefix in [('prodigal','PRODIGAL_RVT'),('sixframe','SIXFRAME_RVT')]:
  rc=run_logged([sys.executable,str(adapter),'--run',str(dep),'--id',cid,'--phase',route],out/(route+'_without_recurrence.log'))
  if rc:raise Hold('historical RT-only should run despite missing recurrence: '+route)
  rr=read_tsv(dep/'scores'/cid/(prefix+'_ONLY.tsv'))
  # Missing recurrence fields differ by design; RT-only rank/evidence must agree.
  signature=lambda rows:[{k:r[k] for k in ['rank','start0','end0','strand','full_evalue','domain_i_evalue','hmm_score']} for r in rows]
  route_checks.append({'method':prefix+'_ONLY','condition':'recurrence_missing',
    'match_original_score_case':signature(rr)==signature(paired_outputs[prefix+'_ONLY'])})
  if (dep/'scores'/cid/(prefix+'_PLUS_FROZEN_REC.tsv')).exists():raise Hold('dependency-missing RT+rec was improperly emitted')
 write_tsv(out/'INDEPENDENT_ROUTE_REGRESSION.tsv',route_checks)
 board=read_tsv(cd/'ARM_M_FULL_FROZEN.tsv');expect=read_tsv(ref/'EXPECTED_ARM_M_CANONICAL_UNION_REGIONS.tsv')
 observed={(int(r['rank']),int(r['start0']),int(r['end0'])) for r in board}
 expected={(int(r['source_rank']),int(r['start0']),int(r['end0'])) for r in expect}
 from historical_compare import compare
 full_equivalence,mismatch=compare(out,ref,f/'vendor')
 for row in route_checks:
  if not row['match_original_score_case']:mismatch.append({'scope':'independent_RT_orchestration','kind':'route_difference','detail':str(row)})
 for r in sorted(expected-observed):mismatch.append({'scope':'ARM_M','kind':'missing_expected_rank_interval','detail':str(r)})
 for r in sorted(observed-expected):mismatch.append({'scope':'ARM_M','kind':'unexpected_rank_interval','detail':str(r)})
 oldrt=read_tsv(ref/'EXPECTED_BASELINE_RT_LOCI.tsv');newrt=read_tsv(cd/'baseline'/cid/'BASELINE_RT_LOCI.tsv')
 key=lambda r:(r['method'],r['strand'],int(r['start0']),int(r['end0']))
 oldmap={key(r):r for r in oldrt};newmap={key(r):r for r in newrt}
 for k in set(oldmap)^set(newmap):mismatch.append({'scope':'RT','kind':'locus_interval_difference','detail':str(k)})
 numerical=[]
 for k in sorted(set(oldmap)&set(newmap)):
  for col in ['full_evalue','domain_i_evalue','hmm_score','domain_score','local_recurrence_hit_n','best_local_recurrence_q']:
   a0=float(oldmap[k][col]);a1=float(newmap[k][col]);same=math.isclose(a0,a1,rel_tol=1e-10,abs_tol=1e-15 if 'q' in col else 0.0)
   numerical.append({'locus':str(k),'column':col,'historical_value':oldmap[k][col],'recomputed_value':newmap[k][col],'match':same})
   if not same:mismatch.append({'scope':'RT','kind':'numeric_difference','detail':str(k)+' '+col})
 write_tsv(out/'HISTORICAL_NUMERIC_COMPARISON.tsv',numerical)
 write_tsv(out/'HISTORICAL_DIFFERENCES.tsv',mismatch,['scope','kind','detail'])
 report={'status':'HISTORICAL_REAL_CPU_REGRESSION_MATCH' if not mismatch else 'HISTORICAL_REAL_CPU_REGRESSION_DIFFERENCES_REQUIRE_REVIEW',
  'sample':'SA1/MW218148.1','whole_native_bp':len(seq),'recomputed_frozen_formal_ACNV_and_recurrence_features':True,
  'legacy_IRR_executed':js(f/'POLICY.json').get('execution_mode','full_legacy')=='full_legacy',
  'execution_mode':js(f/'POLICY.json').get('execution_mode','full_legacy'),
  'full_window_and_board_equivalence':full_equivalence,
  'arm_m_candidate_n':len(board),'expected_arm_m_candidate_n':len(expect),'rt_locus_n':len(newrt),
  'expected_rt_locus_n':len(oldrt),'difference_n':len(mismatch),'numeric_comparison_n':len(numerical),
  'independent_route_checks':len(route_checks),'independent_route_checks_matched':sum(x['match_original_score_case'] for x in route_checks),
  'new_external_benchmark':False,'historically_truth_exposed':True,'scientific_claim':'implementation regression only',
  'source_genome_sequence_sha256':digest(seq),'rvt_hmm_sha256':sha(a.rvt),
  'toolchain':tool,'python_and_package_versions':js(f/'POLICY.json'), 'created_utc':now()}
 write_json(out/'HISTORICAL_REGRESSION_STATUS.json',report)
 write_tsv(out/'HISTORICAL_OUTPUT_SHA256.tsv',hashes(out,[p for p in out.rglob('*') if p.name!='HISTORICAL_OUTPUT_SHA256.tsv']))
 print(json.dumps(report,indent=2));return 0 if not mismatch else 3

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True);ap.add_argument('--cpu',type=int,default=8)
 ap.add_argument('--phase',choices=['all','raw','finish'],default='all');ap.add_argument('--execution-mode',choices=['exact_two_arm','full_legacy'],default='exact_two_arm');ap.add_argument('--rvt');ap.add_argument('--prodigal',default='prodigal');ap.add_argument('--hmmscan',default='hmmscan');ap.add_argument('--hmmpress',default='hmmpress');a=ap.parse_args();out=Path(a.output).resolve()
 try:
  if a.phase!='finish':
   prepare(out,a.cpu,a.execution_mode);arm=out/'scores/HIST_SA1/arm'
   cmd=([sys.executable,str(HERE/'vendor/arm_m/scripts/02_primary_scan.py')] if a.execution_mode=='full_legacy' else
     [sys.executable,str(HERE/'scripts/exact_two_arm_scan.py'),'--vendor',str(out/'frozen/vendor')])
   rc=run_logged(cmd+['--root',str(arm),'--accession','HIST_SA1','--threads',str(a.cpu),'--keep-raw','0'],out/'engine.log')
   if rc:raise Hold('historical full native genome scan failed; inspect engine.log')
  if a.phase!='raw':return finish(a,out)
  print(json.dumps({'status':'HISTORICAL_RAW_SCAN_DONE_BASELINES_PENDING','new_external_validation':False}));return 0
 except (Hold,OSError,ValueError) as ex:print(json.dumps({'status':'HOLD_HISTORICAL_REGRESSION','reason':str(ex),'new_external_validation':False}));return 3
if __name__=='__main__':sys.exit(main())
