#!/usr/bin/env python3
"""Source-locked orchestration; all Stage04A and UG27 inputs are read-only."""
from pathlib import Path,PurePosixPath
import argparse,csv,datetime,gzip,hashlib,json,os,shutil,subprocess,sys,zipfile
import numpy as np
import pandas as pd
from io_utils import *
from scan_groupwise import BOARDS,starts_for,metadata,engine_signature,checkpoint_signature
HERE=Path(__file__).resolve().parents[1]
POL=read_json(HERE/'config/POLICY.json')

def checked_child(root,name):
    p=(root/name).resolve()
    if not p.is_relative_to(root.resolve()):raise RuntimeError('Path outside declared source root: '+name)
    return p

def snapshot(a04):
    p=a04/'config/STAGE04A_POLICY.json';parent=read_json(p)
    z=a04/'provenance'/parent['stage03a_zip']
    check_hash(z,POL['stage03a_zip_sha256'],'frozen engine ZIP')
    roots=[p for p in (a04/'runtime_engine').iterdir() if p.is_dir()]
    if len(roots)!=1:raise RuntimeError('Expected exactly one Stage04A engine directory')
    src=roots[0];dest=HERE/'runtime_engine';dest.mkdir(exist_ok=True)
    records=[]
    with zipfile.ZipFile(z) as zz:
        for name in sorted(zz.namelist()):
            pp=PurePosixPath(name)
            if len(pp.parts)<3 or pp.parts[0]!=src.name:continue
            rel=PurePosixPath(*pp.parts[1:])
            if rel.parts[0] not in ('scripts','config','frozen_sources') or name.endswith('/'):continue
            if '..' in rel.parts or '__pycache__' in rel.parts or str(rel).endswith('.pyc'):continue
            data=zz.read(name);expect=hashlib.sha256(data).hexdigest()
            orig=src/str(rel)
            check_hash(orig,expect,'Stage04A runtime scientific source')
            to=dest/str(rel);to.parent.mkdir(parents=True,exist_ok=True)
            if to.exists():check_hash(to,expect,'existing isolated engine snapshot')
            else:to.write_bytes(data)
            records.append(dict(path=str(rel),sha256=expect))
    if not (dest/'scripts/02_primary_scan.py').exists():raise RuntimeError('No primary entrypoint in frozen ZIP')
    allowed={r['path'] for r in records}
    for p in dest.rglob('*'):
        if p.is_file() and p.parts and '__pycache__' not in p.parts and p.suffix!='.pyc':
            if p.relative_to(dest).as_posix() not in allowed:raise RuntimeError('Extra source file in isolated engine')
    write_tsv(HERE/'outputs/ENGINE_SOURCE_LOCK.tsv',records,['path','sha256'])
    return dest

def run_scan(eng,fasta,acc,out,work,threads,log):
    sig=checkpoint_signature(eng,fasta);g=out/'GROUPWISE_GATE.json'
    if g.exists():
        old=read_json(g)
        if old.get('gate_pass') and old.get('signature')==sig and all((out/n).is_file() and digest(out/n)==old['boards'][n] for n in BOARDS):
            print('[RESUME completed]',acc,flush=True);return old
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[k]='1'
    cmd=[sys.executable,str(HERE/'scripts/scan_groupwise.py'),'--engine',str(eng),'--fasta',str(fasta),
         '--accession',acc,'--output',str(out),'--work',str(work),'--threads',str(threads)]
    log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('a') as f:
        f.write('\nCOMMAND '+json.dumps(cmd)+'\n');f.flush()
        cp=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
    if cp.returncode:raise RuntimeError(f'{acc} replay stopped (rc={cp.returncode}); see {log}')
    return read_json(g)

def byte_compare(ref,new):
    b1=ref.read_bytes();b2=new.read_bytes()
    if ref.name.endswith('.gz'):b1=gzip.decompress(b1);b2=gzip.decompress(b2)
    return dict(reference_sha256=digest(ref),replay_sha256=digest(new),
                serialized_content_byte_identical=(b1==b2),compressed_byte_identical=(ref.read_bytes()==new.read_bytes()))

def verify_validation(a04,eng):
    g=read_json(HERE/'outputs/HISTORICAL_EQUIVALENCE_GATE.json')
    if g.get('gate_pass') is not True:raise RuntimeError('Historical exact-equivalence gate is not PASS')
    if g.get('engine_sha256')!=engine_signature(eng):raise RuntimeError('Engine changed after historical validation')
    if g.get('runtime')!=runtime_info():raise RuntimeError('Runtime changed after validation; repeat validation')
    for name,h in g['implementation_hashes'].items():check_hash(HERE/'scripts'/name,h,'replay implementation')
    for rec in g['comparisons']:
        p=a04/'outputs/genomes'/rec['accession']/rec['board'];check_hash(p,rec['reference_sha256'],'historical board')
        p=HERE/'work/validation'/rec['accession']/'outputs'/rec['board'];check_hash(p,rec['replay_sha256'],'validated replay board')
    return g

def validate(a04,eng,threads):
    check_hash(a04/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv',POL['stage04a_universe_sha256'])
    rows=read_tsv(a04/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv');lookup={r['sequence_id']:r for r in rows}
    comparisons=[]
    for acc in POL['validation_accessions']:
        r=lookup[acc];fasta=checked_child(a04,r['frozen_fasta']);check_hash(fasta,r['fasta_sha256'])
        out=HERE/'work/validation'/acc/'outputs';work=HERE/'work/validation'/acc/'groups'
        print('[HISTORICAL replay]',acc,'log:',HERE/'logs'/f'validate_{acc}.log',flush=True)
        run_scan(eng,fasta,r['analysis_id'],out,work,threads,HERE/'logs'/f'validate_{acc}.log')
        for name in BOARDS:
            ref=a04/'outputs/genomes'/r['analysis_id']/name
            comp=byte_compare(ref,out/name);comp.update(accession=acc,board=name);comparisons.append(comp)
            print(acc,name,comp['serialized_content_byte_identical'],flush=True)
    passed=len(comparisons)==6 and all(r['serialized_content_byte_identical'] for r in comparisons)
    g=dict(stage='Stage04C2-R2C-validation',gate='PASS_HISTORICAL_EXACT_SERIALIZED_BOARDS' if passed else 'HOLD_HISTORICAL_REPLAY_DIFFERENCE',
           gate_pass=passed,comparisons=comparisons,runtime=runtime_info(),engine_sha256=engine_signature(eng),
           implementation_hashes={p.name:digest(p) for p in sorted((HERE/'scripts').glob('*.py'))},
           numerical_tolerance=0.0,ug27_scoring_authorized=passed,truth_read=False,retuning_authorized=False)
    js(HERE/'outputs/HISTORICAL_EQUIVALENCE_GATE.json',g)
    print(json.dumps(g,indent=2))
    if not passed:raise RuntimeError('Historical replay differs; UG27 remains blocked. Do not loosen tolerance.')

def prepare_ug27(a04,ug,eng):
    mf=ug/'outputs/stage04c2_r2_score_freeze/STAGE04C2_R2_CONTIG_SCORING_MANIFEST.tsv'
    rows=read_tsv(mf);want={r['scoring_id']:r for r in POL['contigs']}
    if len(rows)!=15 or len({r['scoring_id'] for r in rows})!=15 or set(want)!={r['scoring_id'] for r in rows}:
        raise RuntimeError('The frozen 15-contig set changed')
    for r in rows:
        w=want[r['scoring_id']]
        if int(r['length_bp'])!=w['length_bp'] or r['sequence_sha256']!=w['sequence_sha256']:raise RuntimeError('Contig content lock differs')
        fp=checked_child(ug,r['frozen_fasta']);check_hash(fp,r['frozen_fasta_file_sha256'],'contig FASTA file')
        header,seq=single_fasta(fp)
        if header!=r['scoring_id'] or len(seq)!=w['length_bp'] or shatext(seq)!=w['sequence_sha256'] or set(seq)-set('ACGT'):
            raise RuntimeError('Contig FASTA semantic lock differs')
    # Verify the full assembly source lock and the decomposition; no locus data.
    smf=ug/'outputs/stage04c2_r2_score_freeze/STAGE04C2_R2_FROZEN_SCORING_SOURCES.tsv'
    check_hash(smf,'48a4d902c0300ae05aa1d5b0ecf1bd76c65df9d0e6fc0c0a500568da111b05e6')
    sources=read_tsv(smf)
    for s in sources:
        p=checked_child(ug,s['native_genome_file']);check_hash(p,s['native_genome_file_sha256'])
        recs=fasta_text(p.read_text());rs=[r for r in rows if r['assembly_id']==s['system_id']]
        if len(recs)!=len(rs):raise RuntimeError('Assembly/contig count mismatch')
        for r in rs:
            seq=recs[r['source_contig_id']]
            if shatext(seq)!=r['sequence_sha256']:raise RuntimeError('Assembly/contig sequence mismatch')
    policy=read_json(eng/'config/ENGINE_POLICY.json');plans=[]
    for r in rows:
        n=int(r['length_bp']);sizes=[len(starts_for(n,L,policy['stride_fraction'])) for L in policy['window_lengths_bp']]
        plans.append(dict(scoring_id=r['scoring_id'],length_bp=n,oriented_windows=2*sum(sizes),
                          largest_stratum_rows=max(sizes),largest_single_memmap_bytes=max(sizes)*5468*4))
    total=sum(r['oriented_windows'] for r in plans)
    if total!=POL['expected_oriented_windows']:raise RuntimeError('Frozen window universe changed')
    write_tsv(HERE/'outputs/RESOURCE_PLAN.tsv',plans,list(plans[0]))
    write_tsv(HERE/'outputs/CONTIG_INPUT_LOCK.tsv',rows,list(rows[0]))
    plan=dict(n_contigs=15,n_assemblies=5,total_bp=sum(r['length_bp'] for r in plans),
              total_oriented_windows=total,cumulative_original_raw_array_bytes=total*(5468+3570)*4,
              original_peak_raw_array_bytes_one_contig=max(r['oriented_windows'] for r in plans)*(5468+3570)*4,
              groupwise_peak_single_memmap_bytes=max(r['largest_single_memmap_bytes'] for r in plans),
              RAM_note='Includes neither robust-score temporary copies nor rank/recurrence DataFrames. Not a RAM upper bound.',
              calibrated_on_complete_strata=True,full_contig_ranking=True)
    js(HERE/'outputs/RESOURCE_PLAN.json',plan);return rows,plan

def ug27(a04,ug,eng,threads):
    historical=verify_validation(a04,eng);rows,plan=prepare_ug27(a04,ug,eng)
    need=2*plan['groupwise_peak_single_memmap_bytes']+2*(1<<30)
    if shutil.disk_usage(HERE).free<need:raise RuntimeError(f'Insufficient working disk; require at least {need} free bytes')
    gpath=HERE/'outputs/UG27_SCORE_FREEZE_GATE.json'
    js(gpath,dict(stage='Stage04C2-R2C',gate='RUNNING_NOT_FROZEN',gate_pass=False,post_score_unblinding_authorized=False))
    completed=[];allrows=[];dual=[]
    for r in sorted(rows,key=lambda x:x['scoring_id']):
        acc=r['scoring_id'];out=HERE/'outputs/contigs'/acc
        print('[UG27]',acc,'log:',HERE/'logs'/f'ug27_{acc}.log',flush=True)
        gg=run_scan(eng,checked_child(ug,r['frozen_fasta']),acc,out,HERE/'work/groups'/acc,threads,HERE/'logs'/f'ug27_{acc}.log')
        completed.append(dict(scoring_id=acc,gate_sha256=digest(out/'GROUPWISE_GATE.json'),**gg['boards']))
        reg=read_tsv(out/BOARDS[2])
        # All raw boards remain unchanged for audit; the public projection only
        # uses the two supported arms, without re-ranking or region re-merging.
        for rr in reg:
            clean={k:v for k,v in rr.items() if k.startswith('acnv_') or k.startswith('recurrence_') or k in ['region_rank','region_id','start0','end0','span_bp','best_q_acnv_family','best_q_recurrence']}
            clean={'assembly_id':r['assembly_id'],'scoring_id':acc,**clean}
            allrows.append(clean)
            fs=set(rr.get('support_families','').split(';'))
            if {'ACNV_FAMILY','GENERIC_RECURRENCE'}<=fs:dual.append(clean)
        js(HERE/'outputs/PROGRESS.json',dict(n_complete=len(completed),n_expected=15,completed=completed,gate_pass=False))
    fields=list(allrows[0]) if allrows else ['assembly_id','scoring_id','region_rank','region_id','start0','end0','span_bp','best_q_acnv_family','best_q_recurrence']
    write_tsv(HERE/'outputs/ALL_REGION_COORDINATES.tsv',allrows,fields)
    write_tsv(HERE/'outputs/TWO_ARM_REGION_VIEW.tsv',dual,fields)
    verify_validation(a04,eng)
    for r in rows:check_hash(checked_child(ug,r['frozen_fasta']),r['frozen_fasta_file_sha256'])
    outman=manifest_tree(HERE/'outputs',HERE/'outputs/SCORE_OUTPUT_MANIFEST.tsv',exclude=['UG27_SCORE_FREEZE_GATE.json'])
    gate=dict(stage='Stage04C2-R2C',gate='PASS_ALL_15_CONTIG_BOARDS_AND_TWO_ARM_VIEW_FROZEN',gate_pass=True,
              n_assemblies=5,n_scoring_contigs=15,n_completed=len(completed),n_all_regions=len(allrows),n_two_arm_regions=len(dual),
              n_oriented_windows=plan['total_oriented_windows'],source_known_score_blind_replay=True,
              truth_used_for_scoring=False,truth_coordinates_loaded=False,retuning_authorized=False,
              score_output_manifest_sha256=outman,historical_equivalence_gate_sha256=digest(HERE/'outputs/HISTORICAL_EQUIVALENCE_GATE.json'),
              engine_signature=engine_signature(eng),post_score_unblinding_authorized=False,
              next_requirement='Collect and archive/hash this score package before a separate locus/array evaluation. No truth read by this runner.')
    js(gpath,gate);print(json.dumps(gate,indent=2))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['validate','plan','score']);ap.add_argument('--stage04a-root',type=Path,required=True);ap.add_argument('--ug27-root',type=Path);ap.add_argument('--threads',type=int,default=8)
    a=ap.parse_args();a04=a.stage04a_root.resolve()
    if a.threads<1:raise RuntimeError('threads must be >=1')
    with lock_dir(HERE/'.run_lock'):
        eng=snapshot(a04)
        if a.mode=='validate':validate(a04,eng,a.threads)
        else:
            if a.ug27_root is None:raise RuntimeError('--ug27-root is required')
            if a.mode=='plan':
                _,plan=prepare_ug27(a04,a.ug27_root.resolve(),eng);print(json.dumps(plan,indent=2))
            else:ug27(a04,a.ug27_root.resolve(),eng,a.threads)
if __name__=='__main__':
    try:main()
    except Exception as e:
        js(HERE/'outputs/LAST_ERROR.json',dict(gate='HOLD_IMPLEMENTATION_OR_INPUT_AUDIT',gate_pass=False,error=str(e),retuning_authorized=False))
        print('[HOLD]',e,file=sys.stderr);raise SystemExit(2)
