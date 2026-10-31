#!/usr/bin/env python3
from __future__ import annotations
import argparse,sys,os,shutil,traceback,platform
from pathlib import Path
from collections import Counter
from common import *
from biology_tools import *
from repeat_evidence import *
from integrate_evidence import integrate
ROOT=Path(__file__).resolve().parents[1]

def preflight(root=ROOT):
    assert_manifest(root,root/'config/INPUT_MANIFEST.tsv')
    c=read_tsv(root/'inputs/candidates.tsv');r=read_tsv(root/'inputs/official_protein_registry.tsv');ss=read_tsv(root/'inputs/genome_sources.tsv');a=fasta_dict(root/'inputs/official_proteins.faa')
    if len(c)!=383 or len({x['sequence_id'] for x in c})!=174 or len(r)!=9576 or len(ss)!=174:raise RuntimeError('Frozen input cardinality mismatch')
    if sha_file(root/'inputs/candidates.tsv')!='73007864fdef7dba662c9b90a02e2b915359c1fbd93a6061edeea35afde07349':raise RuntimeError('Original registry changed')
    if len({x['aa_sha256'] for x in r})!=9330:raise RuntimeError('Unexpected official AA universe')
    if len({x['analysis_id']+'::'+x['region_id'] for x in c})!=383:raise RuntimeError('Duplicate candidate')
    if set(a)!={x['protein_key'] for x in r}:raise RuntimeError('Official peptide key mismatch')
    for x in r:
        seq=a[x['protein_key']][1]
        if sha_text(seq)!=x['aa_sha256'] or len(seq)!=int(x['aa_length']):raise RuntimeError('Official peptide changed')
    lengths={x['accession']:int(x['length_bp']) for x in ss}
    for x in c:
        if not(0<=int(x['start0'])<int(x['end0'])<=lengths[x['sequence_id']]):raise RuntimeError('Candidate coordinate bounds')
    return {'n_candidates':383,'n_accessions':174,'n_official_records':9576,'n_unique_official_AA':len({x['aa_sha256'] for x in r}),'sources_unchanged':True}

def freeze_run(versions):
    ref=read_json(ROOT/'refs/REFERENCE_LOCK.json');assert_manifest(ROOT,ROOT/'refs/REFERENCE_MANIFEST.tsv')
    if ref['method_policy_sha256']!=sha_file(ROOT/'config/METHOD_POLICY.json'):raise RuntimeError('Reference policy mismatch')
    if ref['reference_manifest_sha256']!=sha_file(ROOT/'refs/REFERENCE_MANIFEST.tsv'):raise RuntimeError('Reference manifest mismatch')
    code={str(p.relative_to(ROOT)):sha_file(p) for p in sorted((ROOT/'scripts').glob('*.py'))}
    fp={'input_manifest_sha256':sha_file(ROOT/'config/INPUT_MANIFEST.tsv'),'method_policy_sha256':sha_file(ROOT/'config/METHOD_POLICY.json'),'reference_lock_sha256':sha_file(ROOT/'refs/REFERENCE_LOCK.json'),'script_sha256':code,'tool_versions':versions,'python':platform.python_version()}
    p=ROOT/'outputs/RUN_LOCK.json'
    if p.exists():
        old=read_json(p)
        if old['fingerprint']!=fp:raise RuntimeError('Run inputs/code/tools/references changed. Preserve this run; do not overwrite its lock.')
    else:write_json(p,{'fingerprint':fp,'created_utc':now(),'candidate_retuning_forbidden':True,'UG27_inputs_read':False})
    return sha_file(p)

def real_tool_smoke():
    d=ROOT/'work/tool_smoke';d.mkdir(parents=True,exist_ok=True)
    model=ROOT/'refs/RT_PANEL.hmm';emit=d/'profile_consensus.faa'
    run_cmd(['hmmemit','-c','-o',emit,model],ROOT/'logs/smoke_hmmemit.log')
    targets={x[0]:x[2] for x in fasta(emit)}
    run_cmd(['hmmsearch','--cpu','1','--seed','42','--cut_ga','--noali','--domtblout',d/'consensus.domtbl','-o',d/'consensus.report',model,emit],ROOT/'logs/smoke_hmmsearch.log')
    ds=parse_domtbl(d/'consensus.domtbl',read_tsv(ROOT/'refs/MODEL_CATALOG.tsv'),targets,'GA')
    seen={x['profile_name'] for x in ds}
    required={p['name'] for p in read_json(ROOT/'config/METHOD_POLICY.json')['rt_scan']['profiles']}
    if not required<=seen:raise RuntimeError('Profile-consensus technical smoke failed; do not relax model cutoffs')
    write_json(d/'SMOKE_GATE.json',{'gate_pass':True,'type':'profile consensus technical plumbing test, not biological sensitivity validation','n_profiles':len(required)})

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cpus',type=int,default=4);ap.add_argument('--preflight-only',action='store_true');ap.add_argument('--official-array-preview',type=Path);a=ap.parse_args()
    if a.cpus<1:raise ValueError('cpus must be positive')
    qc=preflight();print('[INPUTS]',json.dumps(qc),flush=True)
    if a.preflight_only:return
    if a.official_array_preview:
        dest=a.official_array_preview;dest.mkdir(parents=True,exist_ok=True)
        rr,ww=build_repeat_evidence(ROOT,None)
        write_tsv(dest/'EXACT_REPEAT_GROUPS.tsv',rr,REPEAT_FIELDS);write_tsv(dest/'RECURRENCE_CODING_CONTEXT.tsv',ww,WINDOW_FIELDS)
        write_json(dest/'PREVIEW_GATE.json',{'gate':'OFFICIAL_ONLY_REPEAT_DESCRIPTORS_COMPUTED_NOT_FINAL_RT_ORF_EVIDENCE','n_candidates':len(ww),'n_repeat_groups':len(rr),'n_candidates_with_repeat_groups':sum(w['n_exact_core_repeat_groups']>0 for w in ww),'method_policy_sha256':sha_file(ROOT/'config/METHOD_POLICY.json'),'HMMER_executed':False,'Prodigal_executed':False,'ART_calls_authorized':False})
        return
    (ROOT/'outputs').mkdir(exist_ok=True);(ROOT/'logs').mkdir(exist_ok=True)
    versions=tool_versions();lock=freeze_run(versions)
    real_tool_smoke()
    all_orfs,near_orfs,aa=run_prodigal_all(ROOT,a.cpus,lock)
    targets,fan=build_targets(ROOT,near_orfs,aa)
    domains=run_hmmer(ROOT,a.cpus,targets,lock)
    rr,ww=build_repeat_evidence(ROOT,all_orfs)
    write_tsv(ROOT/'outputs/EXACT_CORE_REPEAT_GROUPS.tsv',rr,REPEAT_FIELDS);write_tsv(ROOT/'outputs/RECURRENCE_CODING_CONTEXT.tsv',ww,WINDOW_FIELDS)
    gate=integrate(ROOT,fan,domains,rr,ww,all_orfs)
    preflight();assert_manifest(ROOT,ROOT/'refs/REFERENCE_MANIFEST.tsv')
    gate.update({'run_lock_sha256':lock,'sources_unchanged':True,'HMMER_executed':True,'Prodigal_executed':True})
    files=[p for p in (ROOT/'outputs').rglob('*') if p.is_file() and p.name not in ('GATE_STAGE04B_R1D2.json','OUTPUT_MANIFEST.tsv')]
    gate['output_manifest_sha256']=output_manifest(ROOT,files,ROOT/'outputs/OUTPUT_MANIFEST.tsv')
    write_json(ROOT/'outputs/GATE_STAGE04B_R1D2.json',gate);print(json.dumps(gate,indent=2),flush=True)

if __name__=='__main__':
    try: main()
    except Exception as e:
        traceback.print_exc()
        o=ROOT/'outputs';o.mkdir(exist_ok=True)
        # Keep failures visible and never convert missing tools/data into negative biology.
        write_json(o/'GATE_STAGE04B_R1D2.json',{'stage':'Stage04B-R1D2','gate':'HOLD_STAGE04B_R1D2_EXECUTION_OR_INPUT_ERROR','gate_pass':False,'error':f'{type(e).__name__}: {e}','negative_calls_authorized':False,'ART_identity_classification_complete':False,'UG27_touched':False})
        sys.exit(2)
