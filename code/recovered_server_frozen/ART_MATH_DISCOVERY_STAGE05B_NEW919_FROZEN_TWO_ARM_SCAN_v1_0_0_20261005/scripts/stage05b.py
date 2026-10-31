#!/usr/bin/env python3
"""Isolated whole-genome orchestration of the byte-frozen Stage04C1 engine.
No source downloads, protein annotations, biological truth, or detector retuning.
"""
from __future__ import annotations
import argparse, contextlib, csv, datetime, fcntl, gzip, hashlib, importlib.util
import io, json, math, os, platform, shutil, subprocess, sys, tarfile, time
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[1]
ENGINE=BASE/'frozen_engine'
BOARDS=['C1_WINDOW_FAMILY_SCORES.tsv.gz','C1_SELECTED_PROPOSALS.tsv','C1_UNION_REGIONS.tsv']
REG_FIELDS=['region_rank','region_id','start0','end0','span_bp','selected_window_rows','best_q_acnv_family','best_q_recurrence','support_family_n','support_families','sum_two_family_surprisal','best_single_family_q','acnv_best_member','acnv_best_window_start0','acnv_best_window_end0','acnv_best_window_len','recurrence_best_k','recurrence_best_lag_bp','recurrence_best_window_start0','recurrence_best_window_end0','recurrence_best_window_len']

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def json_read(p):return json.loads(Path(p).read_text())
def json_write(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_name(p.name+'.part');q.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');q.replace(p)
def tsv_read(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def tsv_write(p,rows,fields):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_name(p.name+'.part')
    with q.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='raise');w.writeheader();w.writerows(rows)
    q.replace(p)
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def checked(root,relative):
    root=Path(root).resolve();p=(root/relative).resolve()
    require(p.is_relative_to(root),'Path escapes declared source: '+str(relative))
    return p
def file_hash_check(p,h):require(Path(p).is_file() and digest(p)==h,'File missing or SHA256 mismatch: '+str(p))
def runtime():
    return {'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'platform':platform.platform(),'executable':sys.executable}
def pol():return json_read(BASE/'config/STAGE05B_POLICY.json')
def package_check():
    for r in tsv_read(BASE/'PACKAGE_MANIFEST.tsv'):file_hash_check(checked(BASE,r['path']),r['sha256'])
    for r,h in pol()['engine_files'].items():file_hash_check(ENGINE/r,h)
    for r,h in pol()['historical_file_hashes'].items():file_hash_check(BASE/r,h)
def window_count(n):
    return sum(2*((n-L)//(L//4)+1+bool((n-L)%(L//4))) for L in [128,256,512,1024,2048,4096] if n>=L)
def read_fasta(p):
    header=None;parts=[]
    for line in Path(p).read_text(encoding='ascii').splitlines():
        line=line.strip()
        if not line:continue
        if line.startswith('>'):
            require(header is None,'Multi-record FASTA not permitted: '+str(p));header=line[1:]
        else:
            require(header is not None,'Sequence before header: '+str(p));parts.append(line)
    s=''.join(parts)
    require(header and s and not set(s)-set('ACGT'),'Empty/non-ACGT/unnormalized FASTA: '+str(p))
    return header,s

def source_check(source,policy=None,check_all_fastas=True):
    p=policy or pol();source=Path(source).resolve()
    mf=source/'outputs/FROZEN_NEW_GENOME_MANIFEST.tsv'
    file_hash_check(mf,p['source_manifest_sha256'])
    file_hash_check(source/'outputs/SELECTED_SOURCE_ACCESSIONS.tsv',p['selection_sha256'])
    # Freeze gate bytes as well, not just a manually set gate_pass boolean.
    file_hash_check(source/'outputs/SOURCE_COHORT_FREEZE.json',digest(BASE/'config/SOURCE_COHORT_FREEZE.json'))
    file_hash_check(source/'outputs/SELECTION_FREEZE.json',digest(BASE/'config/SELECTION_FREEZE.json'))
    file_hash_check(source/'config/POLICY.json',digest(BASE/'config/STAGE05A_POLICY.json'))
    rr=tsv_read(mf)
    require(len(rr)==p['n_genomes'] and len({r['accession'] for r in rr})==len(rr),'Source cohort count/uniqueness mismatch')
    old={x[k] for x in tsv_read(BASE/'config/EXCLUDE_PRIOR240.tsv') for k in ['sequence_sha256','reverse_complement_sha256']}
    seen=set();audit=[]
    for r in rr:
        a=r['accession'];n=int(r['length_bp'])
        require(r['fasta_path']=='inputs/genomes/'+a+'.frozen.fasta','Unexpected FASTA locator '+a)
        require(window_count(n)==int(r['oriented_windows']),'Window count mismatch '+a)
        canon=min(r['sequence_sha256'],r['reverse_complement_sha256'])
        require(canon not in seen and r['sequence_sha256'] not in old and r['reverse_complement_sha256'] not in old,'Duplicate sequence commitment '+a);seen.add(canon)
        fp=checked(source,r['fasta_path'])
        if check_all_fastas:
            file_hash_check(fp,r['fasta_sha256']);h,s=read_fasta(fp)
            require(h==a and len(s)==n,'Header/length mismatch '+a)
            require(hashlib.sha256(s.encode()).hexdigest()==r['sequence_sha256'],'Sequence hash mismatch '+a)
            rc=s.translate(str.maketrans('ACGT','TGCA'))[::-1]
            require(hashlib.sha256(rc.encode()).hexdigest()==r['reverse_complement_sha256'],'Reverse-complement hash mismatch '+a)
        audit.append({'accession':a,'length_bp':n,'oriented_windows':window_count(n),'fasta_sha256':r['fasta_sha256'],'sequence_sha256':r['sequence_sha256'],'all_byte_checks_pass':check_all_fastas})
    require(sum(int(r['length_bp']) for r in rr)==p['total_bp'],'Total sequence length mismatch')
    require(sum(int(r['oriented_windows']) for r in rr)==p['oriented_windows'],'Total window count mismatch')
    return rr,audit

def build_plan(rr,workers,threads):
    require(workers>0 and threads>0,'Workers and threads must be positive')
    width=5468*4;mat=max(int(r['oriented_windows'])*width for r in rr)
    cpu=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or 1)
    mem=None
    if Path('/proc/meminfo').exists():
        for l in Path('/proc/meminfo').read_text().splitlines():
            if l.startswith('MemAvailable:'):mem=int(l.split()[1])*1024
    disk=shutil.disk_usage(BASE).free;p=pol()
    required_mem=(workers*p['recommended_available_RAM_GiB_per_worker']+p['minimum_extra_available_RAM_GiB'])*2**30
    return {'n_genomes':len(rr),'total_oriented_windows':sum(int(r['oriented_windows']) for r in rr),
            'engine':'BYTE_FROZEN_STAGE04C1_CLEAN_TWO_ARM','max_single_genome_ACNV_memmap_bytes':mat,
            'upper_bound_simultaneous_ACNV_memmap_bytes':min(workers,len(rr))*mat,
            'cumulative_ACNV_array_bytes':sum(int(r['oriented_windows'])*width for r in rr),
            'arrays_released_after_each_successful_genome':True,'memmap_bytes_are_not_RAM_bound':True,
            'genome_workers':workers,'inner_threads':threads,'nominal_feature_worker_processes':workers*threads,
            'available_CPU_affinity':cpu,'available_RAM_bytes_reported_by_OS':mem,'available_disk_bytes':disk,
            'minimum_RAM_bytes_for_this_launch':required_mem,'minimum_disk_bytes_for_this_launch':p['minimum_free_disk_GiB']*2**30,
            'other_running_jobs_not_included_in_CPU_count':True,'no_per_genome_biological_filter':True}

def method_signature():
    files={r['path']:r['sha256'] for r in tsv_read(BASE/'PACKAGE_MANIFEST.tsv') if r['path'].startswith(('scripts/','frozen_engine/','config/'))}
    return hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest()
def work_signature(fasta):return {'method_sha256':method_signature(),'fasta_sha256':digest(fasta),'runtime':runtime()}

def inspect_boards(path,acc,n,expected_fasta=None):
    path=Path(path);g=json_read(path/'C1_SCORE_FREEZE_GATE.json');p=pol()
    require(g.get('gate_pass') is True and g['n_oriented_windows']==window_count(n) and g['length_bp']==n,'Invalid original engine gate '+acc)
    require(g['proposal_families']==p['families'] and g['proposal_q_threshold']==p['proposal_q_threshold'],'Original engine policy mismatch')
    require(g['truth_used'] is False and g['annotation_used'] is False and g['retuning_authorized'] is False,'Unexpected evidence inputs')
    for fn,key in zip(BOARDS,['window_scores_sha256','selected_sha256','union_regions_sha256']):file_hash_check(path/fn,g[key])
    # Validate every window coordinate and finite family score without a new score.
    expected=[(L,st,st+L) for L in [128,256,512,1024,2048,4096] if n>=L for st in starts(n,L)]
    observed=[]
    with gzip.open(path/BOARDS[0],'rt',newline='') as f:
        rd=csv.DictReader(f,delimiter='\t')
        for r in rd:
            require(r['sequence_id']==acc,'Cross-genome window row')
            L,st,en=int(r['window_len']),int(r['start0']),int(r['end0']);observed.append((L,st,en))
            for key in ['q_acnv_family','q_recurrence']:
                v=float(r[key]);require(math.isfinite(v) and 0<v<=1,'Invalid empirical q')
    require(len(observed)==len(expected) and set(observed)==set(expected),'Incomplete/duplicate full-genome windows '+acc)
    regs=tsv_read(path/BOARDS[2]);selected=tsv_read(path/BOARDS[1])
    require(len(regs)==g['n_union_regions'],'Region count mismatch '+acc)
    require([int(r['region_rank']) for r in regs]==list(range(1,len(regs)+1)),'Rank not contiguous '+acc)
    for r in regs:
        require(0<=int(r['start0'])<int(r['end0'])<=n,'Region coordinates outside genome')
        ff=set(r['support_families'].split(';'));require(ff<=set(p['families']) and len(ff)==int(r['support_family_n']),'Unexpected support family')
    for r in selected:require(r['family'] in p['families'] and float(r['family_q'])<=p['proposal_q_threshold'],'Invalid selected proposal')
    return g,regs,selected

def starts(n,L):
    if n<L:return []
    z=list(range(0,n-L+1,L//4))
    if z[-1]!=n-L:z.append(n-L)
    return z

def prepare_job(job):
    for rel in ['config/C1_POLICY.json','frozen_sources/ACNV_SOURCE__legacy_acnv_extractor_frozen.py']:
        src=ENGINE/rel;dst=job/rel;dst.parent.mkdir(parents=True,exist_ok=True)
        if dst.exists():file_hash_check(dst,digest(src))
        else:shutil.copyfile(src,dst)

def execute_genome(acc,fasta,n,job,destination,threads,log):
    sig=work_signature(fasta);done=destination/'COMPLETE.json'
    if destination.exists():
        require(done.exists(),'Output directory without completion marker; preserve for inspection: '+str(destination))
        old=json_read(done);require(old['signature']==sig,'Completed result from different input/method/runtime '+acc)
        for f,h in old['files'].items():file_hash_check(destination/f,h)
        print('[REUSE]',acc,flush=True);return old
    job.mkdir(parents=True,exist_ok=True);prepare_job(job);original=job/'outputs'/acc
    marker=job/'JOB_LOCK.json'
    if marker.exists():require(json_read(marker)['signature']==sig,'Partial job is from a different runtime or source: '+acc)
    else:json_write(marker,{'signature':sig,'accession':acc})
    # Recover a fully completed original-engine run after a wrapper interruption.
    reusable=False
    if (original/'C1_SCORE_FREEZE_GATE.json').exists():
        inspect_boards(original,acc,n);reusable=True
    elapsed=0.0
    if not reusable:
        env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
        for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
        cmd=[sys.executable,str(ENGINE/'scripts/01_scan_two_arm.py'),'--root',str(job),'--control-id',acc,'--accession',acc,'--fasta',str(fasta),'--threads',str(threads),'--keep-raw','0']
        log.parent.mkdir(parents=True,exist_ok=True);t=time.time()
        with log.open('a') as f:
            f.write('\nCOMMAND '+json.dumps(cmd)+'\n');f.flush()
            cp=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
        elapsed=time.time()-t
        require(cp.returncode==0,f'{acc}: original frozen engine failed rc={cp.returncode}; see {log}')
    file_hash_check(fasta,sig['fasta_sha256'])
    g,regs,selected=inspect_boards(original,acc,n)
    info={'accession':acc,'signature':sig,'length_bp':n,'n_oriented_windows':window_count(n),'n_union_regions':len(regs),
          'n_dual_arm_regions':sum(int(r['support_family_n'])==2 for r in regs),'elapsed_seconds_this_attempt':elapsed,
          'feature_seconds_original_engine':g['acnv_seconds'],'threads_this_attempt':threads,
          'files':{f:digest(original/f) for f in BOARDS+['C1_SCORE_FREEZE_GATE.json']},'gate_pass':True,
          'biological_annotation_authorized':False,'all_union_regions_retained':True,'mathematical_rank_is_not_ART_identity':True}
    destination.parent.mkdir(parents=True,exist_ok=True)
    tmp=destination.with_name(destination.name+'.partial_commit')
    require(not tmp.exists(),'Partial commit exists: preserve and inspect '+str(tmp))
    shutil.copytree(original,tmp);json_write(tmp/'COMPLETE.json',info);tmp.rename(destination)
    # The source engine normally removes this. Recovery may have interrupted
    # that last housekeeping line; only this package's own numeric work file.
    raw=job/'work'/acc/'acnv_5468.npy'
    if raw.exists():raw.unlink()
    return info

def engine_import():
    sys.path.insert(0,str(ENGINE/'scripts'))
    sp=importlib.util.spec_from_file_location('frozen_two_arm_reference',ENGINE/'scripts/01_scan_two_arm.py')
    mod=importlib.util.module_from_spec(sp);sys.modules[sp.name]=mod;sp.loader.exec_module(mod)
    return mod

def validate(threads):
    gate=BASE/'outputs/HISTORICAL_EQUIVALENCE_GATE.json';signature={'method_sha256':method_signature(),'runtime':runtime()}
    if gate.exists():
        old=json_read(gate)
        require(old.get('gate_pass') is True,'Previous validation was HOLD. Preserve it; do not bypass.')
        require(old['signature']==signature,'Validation/runtime changed; use a separate package copy.')
        for row in old['comparisons']:
            for f,h in row['actual_files'].items():file_hash_check(BASE/f,h)
        print('[REUSE] historical validation exact PASS',flush=True);return old
    comparisons=[];mod=engine_import();ep=json_read(ENGINE/'config/C1_POLICY.json')
    for a in pol()['validation_accessions']:
        fasta=BASE/'validation_inputs'/f'{a}.frozen.fasta';h,s=read_fasta(fasta)
        require(h.split()[0]==a,'Historical identifier mismatch')
        dst=BASE/'work/validation'/a/'verified';job=BASE/'work/validation'/a/'job'
        print('[HISTORICAL]',a,'log=',BASE/'logs'/('historical_'+a+'.log'),flush=True)
        execute_genome(a,fasta,len(s),job,dst,threads,BASE/'logs'/('historical_'+a+'.log'))
        expected_path=BASE/'validation_expected'/f'{a}.STAGE04A_WINDOW_FAMILY_SCORES.tsv.gz'
        text=gzip.decompress(expected_path.read_bytes()).decode()
        rows=list(csv.reader(io.StringIO(text),delimiter='\t'));head=rows[0]
        with gzip.open(dst/BOARDS[0],'rt',newline='') as f:new=list(csv.reader(f,delimiter='\t'))
        ix=[head.index(c) for c in new[0]]
        # All common fields serialized in the original old genome score board,
        # INCLUDING every ACNV member score, strand/tie and recurrence field.
        expected=[[r[i] for i in ix] for r in rows]
        common_exact=(new==expected)
        d=pd.read_csv(expected_path,sep='\t',float_precision='round_trip')
        chunks=[]
        for fam,qc in [('ACNV_FAMILY','q_acnv_family'),('GENERIC_RECURRENCE','q_recurrence')]:
            z=d.loc[d[qc]<=ep['proposal_q_threshold'],['sequence_id','window_len','start0','end0',qc]].copy();z['family']=fam;z=z.rename(columns={qc:'family_q'});chunks.append(z)
        selected=pd.concat(chunks,ignore_index=True).sort_values(['start0','end0','family'])
        regions=mod.merge_two(selected,d)
        selected_exact=selected.to_csv(sep='\t',index=False).encode()==(dst/BOARDS[1]).read_bytes()
        regions_exact=regions.to_csv(sep='\t',index=False).encode()==(dst/BOARDS[2]).read_bytes()
        z={'accession':a,'common_window_fields_exact':common_exact,'clean_two_arm_proposals_serialized_exact':selected_exact,
           'clean_two_arm_regions_serialized_exact':regions_exact,'expected_regions_derived_from_old_frozen_scores_not_old_pilot_ranks':True,
           'actual_files':{str((dst/f).relative_to(BASE)):digest(dst/f) for f in BOARDS},'n_genomic_windows':len(new)-1}
        comparisons.append(z);print('[HISTORICAL RESULT]',json.dumps(z),flush=True)
    passed=all(r['common_window_fields_exact'] and r['clean_two_arm_proposals_serialized_exact'] and r['clean_two_arm_regions_serialized_exact'] for r in comparisons)
    out={'stage':'Stage05B-validation','gate':'PASS_HISTORICAL_COMMON_SCORES_AND_CLEAN_TWO_ARM_BOARDS_EXACT' if passed else 'HOLD_HISTORICAL_EQUIVALENCE','gate_pass':passed,'signature':signature,'comparisons':comparisons,
         'new919_scored_during_validation':False,'old240_or_383_changed':False,'UG27_read_or_written':False,'numeric_tolerance_for_PASS':0}
    json_write(gate,out);require(passed,'Historical equivalence failed; no new919 scoring authorized')
    return out

def preflight(source,workers,threads):
    package_check();rr,audit=source_check(source);plan=build_plan(rr,workers,threads)
    json_write(BASE/'outputs/RESOURCE_PLAN.json',plan)
    tsv_write(BASE/'outputs/SOURCE919_RECHECK.tsv',audit,list(audit[0]))
    require(plan['nominal_feature_worker_processes']+workers<=plan['available_CPU_affinity'],'Too many local workers for CPU affinity; lower workers/threads without changing sources')
    require(plan['available_disk_bytes']>=plan['minimum_disk_bytes_for_this_launch'],'Insufficient disk for launch; see RESOURCE_PLAN.json')
    if plan['available_RAM_bytes_reported_by_OS'] is not None:
        require(plan['available_RAM_bytes_reported_by_OS']>=plan['minimum_RAM_bytes_for_this_launch'],'Insufficient reported available RAM for launch; reduce genome-workers, not cohort')
    lock={'stage':'Stage05B','source_root':str(Path(source).resolve()),'source_manifest_sha256':pol()['source_manifest_sha256'],'selection_sha256':pol()['selection_sha256'],'n_genomes':len(rr),'method_sha256':method_signature(),'runtime':runtime()}
    dest=BASE/'outputs/RUN_LOCK.json'
    if dest.exists():require(json_read(dest)==lock,'Frozen run lock differs; do not mix runtimes, paths or methods')
    else:json_write(dest,lock)
    json_write(BASE/'outputs/PREFLIGHT_GATE.json',{'gate':'PASS_ALL919_SOURCE_BYTES_AND_FROZEN_ENGINE_VERIFIED','gate_pass':True,'n_frozen_fastas_rehashed':len(rr),'n_genomes':len(rr),'annotation_used':False,'UG27_touched':False,'plan':plan})
    print('[PREFLIGHT PASS]',len(rr),'source FASTAs verified; max ACNV file bytes=',plan['max_single_genome_ACNV_memmap_bytes'],flush=True)
    return rr

def run_all(source,workers,threads):
    rr=preflight(source,workers,threads);validate(threads)
    queue=sorted(rr,key=lambda r:(-int(r['length_bp']),r['accession']))
    tsv_write(BASE/'outputs/FROZEN_EXECUTION_ORDER.tsv',[{'execution_order':i+1,'accession':r['accession'],'length_bp':r['length_bp'],'oriented_windows':r['oriented_windows']} for i,r in enumerate(queue)],['execution_order','accession','length_bp','oriented_windows'])
    def one(r):
        a=r['accession'];fp=checked(source,r['fasta_path']);file_hash_check(fp,r['fasta_sha256'])
        return execute_genome(a,fp,int(r['length_bp']),BASE/'work/jobs'/a,BASE/'outputs/genomes'/a,threads,BASE/'logs'/('scan_'+a+'.log'))
    start=time.time();state={'phase':'SCANNING_ALL919','completed_verified_in_this_invocation':0,'cohort_total':len(rr),'annotation_authorized':False};json_write(BASE/'outputs/RUN_STATUS.json',state)
    itr=iter(queue);futures={};errors=[];n_done=0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in range(workers):
            r=next(itr,None)
            if r is not None:futures[pool.submit(one,r)]=r
        while futures:
            ready,_=wait(list(futures),return_when=FIRST_COMPLETED)
            for fu in ready:
                r=futures.pop(fu)
                try:
                    z=fu.result();n_done+=1;print(f"[{n_done}/{len(rr)}] SCORE_FROZEN {r['accession']} regions={z['n_union_regions']} dual={z['n_dual_arm_regions']}",flush=True)
                except Exception as e:errors.append({'accession':r['accession'],'error':str(e)});print('[HOLD]',r['accession'],str(e),flush=True)
                state.update(completed_verified_in_this_invocation=n_done,errors=errors,elapsed_seconds=time.time()-start);json_write(BASE/'outputs/RUN_STATUS.json',state)
                if not errors:
                    nr=next(itr,None)
                    if nr is not None:futures[pool.submit(one,nr)]=nr
    if errors:
        json_write(BASE/'outputs/STAGE05B_FINAL_GATE.json',{'gate':'HOLD_STAGE05B_INCOMPLETE_SCORING','gate_pass':False,'errors':errors,'n_expected':len(rr),'annotation_authorized':False})
        raise RuntimeError('Scoring stopped safely. Completed genomes retained; fix the reported cause and resume the same command.')
    finalize(source)

def finalize(source):
    package_check();rr,_=source_check(source);lock=json_read(BASE/'outputs/RUN_LOCK.json')
    require(lock['method_sha256']==method_signature() and lock['runtime']==runtime(),'Method or runtime changed')
    require(json_read(BASE/'outputs/HISTORICAL_EQUIVALENCE_GATE.json')['gate_pass'] is True,'Historical validation missing')
    allrows=[];burden=[];outfiles=[]
    for r in rr:
        a=r['accession'];dst=BASE/'outputs/genomes'/a;info=json_read(dst/'COMPLETE.json')
        require(info['signature']==work_signature(checked(source,r['fasta_path'])),'Per-genome signature mismatch '+a)
        for f,h in info['files'].items():file_hash_check(dst/f,h);outfiles.append({'path':str((dst/f).relative_to(BASE)),'sha256':h,'size':(dst/f).stat().st_size})
        regs=tsv_read(dst/BOARDS[2]);nd=0
        for z in regs:
            dual=set(z['support_families'].split(';'))==set(pol()['families']);nd+=int(dual)
            allrows.append(dict(sequence_id=a,analysis_id=a,cohort='STAGE05_NEW919',source_bin=r['bin'],length_bp=r['length_bp'],**z,dual_arm_candidate=int(dual)))
        burden.append({'sequence_id':a,'source_bin':r['bin'],'length_bp':r['length_bp'],'n_union_regions':len(regs),'n_dual_arm_regions':nd,'n_oriented_windows':r['oriented_windows'],'all_sources_retained':True})
    cols=['sequence_id','analysis_id','cohort','source_bin','length_bp']+REG_FIELDS+['dual_arm_candidate']
    allp=BASE/'outputs/STAGE05B_ALL_TWO_ARM_UNION_REGIONS.tsv';dp=BASE/'outputs/STAGE05B_DUAL_ARM_CANDIDATES.tsv';bp=BASE/'outputs/STAGE05B_PER_GENOME_BURDEN.tsv'
    tsv_write(allp,allrows,cols);tsv_write(dp,[r for r in allrows if r['dual_arm_candidate']],cols);tsv_write(bp,burden,list(burden[0]))
    for p in [allp,dp,bp]:outfiles.append({'path':str(p.relative_to(BASE)),'sha256':digest(p),'size':p.stat().st_size})
    tsv_write(BASE/'outputs/SCORE_FILE_MANIFEST.tsv',outfiles,['path','sha256','size'])
    gate={'stage':'Stage05B','gate':'PASS_ALL919_NEW_GENOME_TWO_ARM_SCORE_BOARDS_FROZEN','gate_pass':True,
          'n_expected_genomes':len(rr),'n_success_genomes':len(rr),'n_failed_genomes':0,'total_bp':pol()['total_bp'],'oriented_windows':pol()['oriented_windows'],
          'n_all_union_regions':len(allrows),'n_dual_arm_candidates':sum(r['dual_arm_candidate'] for r in allrows),'n_genomes_with_dual_arm':sum(r['n_dual_arm_regions']>0 for r in burden),
          'source_manifest_sha256':pol()['source_manifest_sha256'],'all_regions_sha256':digest(allp),'dual_candidates_sha256':digest(dp),'score_manifest_sha256':digest(BASE/'outputs/SCORE_FILE_MANIFEST.tsv'),
          'run_lock_sha256':digest(BASE/'outputs/RUN_LOCK.json'),'all_919_genomes_retained':True,'no_topk_filter':True,
          'annotation_used_in_scoring':False,'new_ART_systems_confirmed':0,'new_to_project_not_publication_prospective':True,
          'near_homology_independence_established':False,'old240_modified':False,'old383_modified':False,'UG27_touched':False,
          'retuning_authorized':False,'post_score_architecture_and_novelty_stage_ready':True,'archive_results_and_sha_before_biology':True,
          'engine_lineage':'Existing Stage04C1 clean two-arm byte-frozen source; legacy pilot region identities are not redefined.'}
    json_write(BASE/'outputs/STAGE05B_FINAL_GATE.json',gate);json_write(BASE/'outputs/RUN_STATUS.json',{'phase':'ALL919_SCORE_FREEZE_COMPLETE','n_completed':len(rr),'n_total':len(rr),'biology_performed':False})
    print(json.dumps(gate,indent=2),flush=True)

def status():
    rr=tsv_read(BASE/'config/FROZEN_NEW_GENOME_MANIFEST.tsv');done=list((BASE/'outputs/genomes').glob('*/COMPLETE.json'))
    print('Completed markers:',len(done),'/',len(rr),'(status view; finalization checks file hashes)')
    for f in ['RUN_STATUS.json','STAGE05B_FINAL_GATE.json']:
        p=BASE/'outputs'/f
        if p.exists():print(p.name);print(p.read_text())

def collect(with_windows=False):
    # Refuse collection while the run owns the write lock.
    package_check();stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S');items=[];omitted=[]
    for folder in ['outputs','logs','config','docs','evidence','scripts','tests','frozen_engine','validation_inputs','validation_expected']:
        for p in (BASE/folder).rglob('*'):
            if not p.is_file() or '__pycache__' in p.parts or p.suffix=='.pyc' or p.name.endswith('.part'):continue
            if p.name==BOARDS[0] and folder=='outputs' and not with_windows:
                omitted.append({'path':str(p.relative_to(BASE)),'sha256':digest(p),'size':p.stat().st_size});continue
            items.append(p)
    for p in BASE.iterdir():
        if p.is_file() and p.suffix in ['.sh','.md','.tsv']:items.append(p)
    om=BASE/'evidence/RESULT_COLLECTION_SCOPE.json';json_write(om,{'window_score_bulk_included':with_windows,'omitted_window_boards':omitted,'raw_FASTAs_919_included':False,'raw_arrays_included':False,'partial_results_allowed':True});items.append(om)
    dest=BASE/('ART_STAGE05B_NEW919_SCORE_RESULTS_'+stamp+'.tar.gz');require(not dest.exists(),'Collection name exists; retry after one second')
    with tarfile.open(dest,'w:gz') as t:
        for p in sorted(set(items)):t.add(p,arcname=str(p.relative_to(BASE)),recursive=False)
    Path(str(dest)+'.sha256').write_text(digest(dest)+'  '+dest.name+'\n');print(dest);print(digest(dest))

@contextlib.contextmanager
def exclusive():
    p=BASE/'.stage05b.flock';f=p.open('a+')
    try:fcntl.flock(f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:f.close();raise RuntimeError('Another Stage05B writer/collector is active; do not launch duplicate jobs')
    try:
        f.seek(0);f.truncate();f.write(str(os.getpid())+'\n');f.flush();yield
    finally:fcntl.flock(f.fileno(),fcntl.LOCK_UN);f.close()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['preflight','validate','run','finalize','status','collect']);ap.add_argument('--source',type=Path);ap.add_argument('--genome-workers',type=int,default=2);ap.add_argument('--threads',type=int,default=8);ap.add_argument('--with-window-scores',action='store_true');a=ap.parse_args()
    if a.action=='status':status();return
    (BASE/'outputs').mkdir(exist_ok=True)
    with exclusive():
        package_check()
        if a.action in ['preflight','run','finalize']:require(a.source is not None,'--source is required')
        if a.action=='preflight':preflight(a.source,a.genome_workers,a.threads)
        elif a.action=='validate':validate(a.threads)
        elif a.action=='run':run_all(a.source,a.genome_workers,a.threads)
        elif a.action=='finalize':finalize(a.source)
        elif a.action=='collect':collect(a.with_window_scores)
if __name__=='__main__':
    try:main()
    except Exception as e:
        print('[FAIL-CLOSED]',type(e).__name__,str(e),file=sys.stderr);sys.exit(2)
