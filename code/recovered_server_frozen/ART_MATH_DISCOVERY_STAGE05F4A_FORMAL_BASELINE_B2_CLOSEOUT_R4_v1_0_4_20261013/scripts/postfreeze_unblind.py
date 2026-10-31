#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
from common import *

def ov(a,b,c,d): return a<d and b>c

def truth_rank(rows,t0,t1,start='eval_start0',end='eval_end0'):
    for i,r in enumerate(rows,1):
        if ov(int(r[start]),int(r[end]),t0,t1): return i
    return None

def rt_orders(rows,method,n_bp,radius):
    x=[r.copy() for r in rows if r['method']==method]
    for r in x:
        s=max(0,int(r['start0'])-radius); e=min(n_bp,int(r['end0'])+radius)
        r['eval_start0']=s; r['eval_end0']=e
    rt=sorted(x,key=lambda r:(float(r['full_evalue']),float(r['domain_i_evalue']),-float(r['hmm_score']),int(r['start0'])))
    rr=sorted(x,key=lambda r:(0 if int(r['local_recurrence_hit_n'])>0 else 1,float(r['best_local_recurrence_q']),float(r['full_evalue']),float(r['domain_i_evalue']),-float(r['hmm_score']),int(r['start0'])))
    return rt,rr

def arm_rows(rows,start='start0',end='end0'):
    z=[]
    for r in rows:
        q=r.copy(); q['eval_start0']=q[start]; q['eval_end0']=q[end]; z.append(q)
    return z

def add(out,c,method,arr,t0,t1,cdef,burden_per_mb='',exposure_fraction='',truth_dual_arm=''):
    r=truth_rank(arr,t0,t1); n=len(arr)
    out.append({
        'case_id':c['case_id'],'group':c['group'],'role':c['role'],'truth_mode':c.get('truth_mode',''),
        'truth_independence_tier':c.get('truth_independence_tier',''),'method':method,'candidate_definition':cdef,
        'candidate_n':n,'truth_rank':r or '','recovered':int(r is not None),'hit1':int(r==1),
        'hit3':int(r is not None and r<=3),'hit5':int(r is not None and r<=5),'hit10':int(r is not None and r<=10),
        'reciprocal_rank':1.0/r if r else 0.0,'normalized_rank_percentile':100.0*r/n if r and n else '',
        'candidate_burden_per_mb':burden_per_mb,'genome_exposure_fraction':exposure_fraction,
        'truth_dual_arm':truth_dual_arm,
    })

def summarize(z,m):
    if not z: return None
    n=len(z); rec=sum(int(r['recovered']) for r in z)
    def mean(key): return sum(float(r[key]) for r in z)/n
    ranks=[int(r['truth_rank']) for r in z if str(r['truth_rank'])!='']
    return {'method':m,'n':n,'recovered_n':rec,'recovered_fraction':rec/n,'hit1_fraction':mean('hit1'),'hit3_fraction':mean('hit3'),'hit5_fraction':mean('hit5'),'hit10_fraction':mean('hit10'),'mrr':mean('reciprocal_rank'),'mean_normalized_rank_percentile':sum(float(r['normalized_rank_percentile']) for r in z if str(r['normalized_rank_percentile'])!='')/len(ranks) if ranks else '','mean_candidate_burden_per_mb':mean('candidate_burden_per_mb'),'mean_genome_exposure_fraction':mean('genome_exposure_fraction')}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--truth-tsv',required=True); a=ap.parse_args()
    here=Path(__file__).resolve().parents[1]
    contract=json.loads((here/'config/BENCHMARK_CONTRACT.json').read_text()); radius=int(contract['primary_rt_neighborhood_radius_bp'])
    b1=json.loads((here/'outputs/score_freeze/BASELINE_B1_GATE.json').read_text())
    if not b1['gate_pass']: raise SystemExit('B1 not frozen')
    if b1.get('truth_opened_for_scoring') is not False: raise SystemExit('B1 truth chronology invalid')
    cases=read_tsv(Path(a.truth_tsv))
    unresolved=[c['case_id'] for c in cases if c.get('primary_eval')=='1' and (c.get('truth_start0') in {'','.','NA'} or c.get('truth_end0') in {'','.','NA'})]
    if unresolved:
        gate={'stage':'Stage05F4A-B2-R4-postfreeze-evaluation','version':'v1.0.4','gate':'HOLD_STAGE05F4A_B2_R4_UNRESOLVED_PROJECT_TRUTH','gate_pass':False,'unresolved_truth_cases':unresolved,'truth_used_before_score_freeze':False,'retuning_authorized':False}
        safe_json(here/'outputs/BASELINE_B2_GATE.json',gate); print(json.dumps(gate,indent=2)); raise SystemExit(2)
    seq={r['case_id']:r for r in read_tsv(here/'frozen/BASELINE_SEQUENCE_MANIFEST.tsv')}
    summ={r['case_id']:r for r in read_tsv(here/'outputs/score_freeze/BASELINE_SCORE_SUMMARY.tsv')}
    out=[]
    for c in cases:
        if c.get('primary_eval')!='1': continue
        t0=int(c['truth_start0']); t1=int(c['truth_end0']); cid=c['case_id']; n_bp=int(seq[cid]['length_bp']); mb=n_bp/1e6
        props=read_tsv(here/'frozen/canonical'/cid/'CANONICAL_PROPOSALS.tsv')
        union_raw=sorted(read_tsv(here/'frozen/canonical'/cid/'CANONICAL_UNION_REGIONS.tsv'),key=lambda r:(int(r['source_rank']),int(r['source_row'])))
        truth_dual=int(any(ov(int(r['start0']),int(r['end0']),t0,t1) and int(r.get('dual_arm','0'))==1 for r in union_raw))
        unions=arm_rows(union_raw)
        ac=arm_rows(sorted([r for r in props if 'ACNV' in r['family']],key=lambda r:(float(r['qvalue']),int(r['source_row']))))
        rr=arm_rows(sorted([r for r in props if 'RECURRENCE' in r['family']],key=lambda r:(float(r['qvalue']),int(r['source_row']))))
        add(out,c,'ARM_M_FULL_FROZEN',unions,t0,t1,'frozen union region',len(unions)/mb if mb else '',1.0,truth_dual)
        add(out,c,'ACNV_ONLY_FROZEN',ac,t0,t1,'frozen ACNV proposal window',len(ac)/mb if mb else '',1.0,'')
        add(out,c,'RECURRENCE_ONLY_FROZEN',rr,t0,t1,'frozen recurrence proposal window',len(rr)/mb if mb else '',1.0,'')
        rows=read_tsv(here/'outputs/score_freeze'/cid/'BASELINE_RT_LOCI.tsv')
        for method,prefix in [('PRODIGAL_RVT','prodigal_rvt'),('SIXFRAME_RVT','sixframe_rvt')]:
            rt,rrec=rt_orders(rows,method,n_bp,radius)
            burden=float(summ[cid][prefix+'_locus_per_mb']); exposure=float(summ[cid][prefix+'_exposure_fraction'])
            add(out,c,method+'_ONLY',rt,t0,t1,f'RT-centered +/-{radius} bp neighborhood ranked by frozen HMM evidence',burden,exposure,'')
            add(out,c,method+'_PLUS_FROZEN_REC',rrec,t0,t1,f'same RT-centered +/-{radius} bp neighborhood ranked by frozen recurrence then HMM evidence',burden,exposure,'')
    write_tsv(here/'outputs/BASELINE_POSTFREEZE_EVALUATION.tsv',out,list(out[0]))
    methods=[]
    for r in out:
        if r['method'] not in methods: methods.append(r['method'])
    primary=[]; exact=[]; strict=[]; reconstructed=[]
    for m in methods:
        z=[r for r in out if r['method']==m and r['group']=='ART_POSITIVE' and r['role']=='PRIMARY_EXTERNAL']
        s=summarize(z,m)
        if s: primary.append(s)
        ze=[r for r in z if r['truth_mode'] in {'STRICT_REPEAT_CORE','PROJECT_POSTFREEZE_EXACT_REPEAT_ENVELOPE'}]
        s=summarize(ze,m)
        if s: exact.append(s)
        zs=[r for r in z if r['truth_mode']=='STRICT_REPEAT_CORE']
        s=summarize(zs,m)
        if s: strict.append(s)
        zr=[r for r in z if r['truth_mode']=='PROJECT_POSTFREEZE_EXACT_REPEAT_ENVELOPE']
        s=summarize(zr,m)
        if s: reconstructed.append(s)
    if primary: write_tsv(here/'outputs/BASELINE_PRIMARY_SUMMARY.tsv',primary,list(primary[0]))
    if exact: write_tsv(here/'outputs/BASELINE_EXACT_ARRAY_SUBSET_SUMMARY.tsv',exact,list(exact[0]))
    if strict: write_tsv(here/'outputs/BASELINE_STRICT_EXTERNAL_ARRAY_SUBSET_SUMMARY.tsv',strict,list(strict[0]))
    if reconstructed: write_tsv(here/'outputs/BASELINE_PROJECT_RECONSTRUCTED_ARRAY_SUBSET_SUMMARY.tsv',reconstructed,list(reconstructed[0]))
    ctrl=[]
    arm=[r for r in out if r['group']=='NONARRAY_RT_CONTROL' and r['method']=='ARM_M_FULL_FROZEN']
    for r in arm: ctrl.append({'case_id':r['case_id'],'truth_rank':r['truth_rank'],'candidate_n':r['candidate_n'],'normalized_rank_percentile':r['normalized_rank_percentile'],'truth_dual_arm':r['truth_dual_arm']})
    if ctrl: write_tsv(here/'outputs/BASELINE_CONTROL_DUAL_SIGNAL.tsv',ctrl,list(ctrl[0]))
    dual=sum(int(r['truth_dual_arm']) for r in arm) if arm else 0
    gate={'stage':'Stage05F4A-B2-R4-postfreeze-evaluation','version':'v1.0.4','gate':'PASS_STAGE05F4A_B2_R4_FORMAL_BASELINE_POSTFREEZE_EVALUATION_COMPLETE','gate_pass':True,'truth_used_before_score_freeze':False,'retuning_authorized':False,'protein_first_candidate_definition':f'prespecified +/-{radius} bp RT-centered neighborhood; ranking frozen in B1','pre_registered_primary_positive_n':len([c for c in cases if c.get('group')=='ART_POSITIVE' and c.get('role')=='PRIMARY_EXTERNAL']),'exact_array_descriptive_subset_n':len([c for c in cases if c.get('group')=='ART_POSITIVE' and c.get('role')=='PRIMARY_EXTERNAL' and c.get('truth_mode') in {'STRICT_REPEAT_CORE','PROJECT_POSTFREEZE_EXACT_REPEAT_ENVELOPE'}]),'strict_external_array_subset_n':len([c for c in cases if c.get('group')=='ART_POSITIVE' and c.get('role')=='PRIMARY_EXTERNAL' and c.get('truth_mode')=='STRICT_REPEAT_CORE']),'project_reconstructed_array_subset_n':len([c for c in cases if c.get('group')=='ART_POSITIVE' and c.get('role')=='PRIMARY_EXTERNAL' and c.get('truth_mode')=='PROJECT_POSTFREEZE_EXACT_REPEAT_ENVELOPE']),'nonarray_rt_control_n':len(arm),'nonarray_rt_control_dual_signal_n':dual,'nonarray_rt_control_dual_signal_fraction':dual/len(arm) if arm else None,'note':'B0/B1 scores inherited unchanged from archived R2. The 8-case primary panel is preserved; the 7-case exact-array subset is descriptive sensitivity. PB50 remains the prelocked RT-proxy challenge.'}
    safe_json(here/'outputs/BASELINE_B2_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
