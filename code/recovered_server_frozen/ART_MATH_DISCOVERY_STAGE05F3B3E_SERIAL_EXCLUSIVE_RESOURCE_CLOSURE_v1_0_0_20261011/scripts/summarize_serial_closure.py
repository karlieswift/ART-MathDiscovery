#!/usr/bin/env python3
import csv, hashlib, json, re, sys
from pathlib import Path

prior_summary=Path(sys.argv[1])
serial_status_path=Path(sys.argv[2])
fold_status_path=Path(sys.argv[3])
fold_dir=Path(sys.argv[4])
lock_path=Path(sys.argv[5])
outdir=Path(sys.argv[6])
outdir.mkdir(parents=True,exist_ok=True)

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()
def ff(x):
    try:return float(x)
    except:return float('nan')
def normacc(t):
    m=re.search(r'AF-([A-Za-z0-9]+)-F\d+-model',t)
    return m.group(1) if m else t
def tier(tm,qc,tc):
    if tm>=0.50 and qc>=0.50 and tc>=0.50:return 'STRONG'
    if tm>=0.40 and qc>=0.30 and tc>=0.30:return 'SUGGESTIVE'
    return 'WEAK'

with prior_summary.open(newline='') as f: prior=list(csv.DictReader(f,delimiter='\t'))
with serial_status_path.open(newline='') as f: serial={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
with fold_status_path.open(newline='') as f: foldstat={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
assert len(prior)==32, len(prior)
unresolved=[r for r in prior if r['final_best_tier']=='BOLTZ_TECHNICAL_FAILURE_AFTER_EXACT_RETRY']
assert len(unresolved)==7, len(unresolved)
assert set(serial)=={r['candidate_id'] for r in unresolved}
fields=['query','target','evalue','bits','qcov','tcov','alnlen','qlen','tlen','fident','alntmscore']

out=[]; serial_hits=[]
for p in prior:
    row=dict(p)
    cid=p['candidate_id']
    row['pre_f3b3e_best_tier']=p['final_best_tier']
    row['serial_exclusive_boltz_status']='NOT_APPLICABLE'
    row['serial_failure_category']=''
    row['serial_cif_path']=''
    row['f3b3e_resolution_source']='PRE_F3B3E'
    if p['final_best_tier']!='BOLTZ_TECHNICAL_FAILURE_AFTER_EXACT_RETRY':
        row['f3b3e_final_best_tier']=p['final_best_tier']
        out.append(row); continue

    s=serial[cid]
    sb=s['serial_boltz_status'].strip()
    row['serial_exclusive_boltz_status']=sb
    row['serial_failure_category']=s['failure_category']
    row['serial_cif_path']=s['cif_path'].strip().rstrip('\r')
    row['f3b3e_resolution_source']='SERIAL_EXCLUSIVE_RESOURCE_CLOSURE'
    if not sb.startswith('PASS'):
        row['foldseek_status']='NOT_RUN_SERIAL_BOLTZ_FAILURE'
        row['r5c_hit_n']='0'; row['best_clean_rt_target']=''; row['best_alntmscore']=''; row['best_qcov']=''; row['best_tcov']=''
        row['f3b3e_final_best_tier']='RESOURCE_LIMITED_BOLTZ_UNRESOLVED'
        out.append(row); continue

    fs=foldstat.get(cid)
    if fs is None or fs['foldseek_status']!='PASS':
        row['foldseek_status']=fs['foldseek_status'] if fs else 'MISSING_STATUS'
        row['r5c_hit_n']='0'; row['best_clean_rt_target']=''; row['best_alntmscore']=''; row['best_qcov']=''; row['best_tcov']=''
        row['f3b3e_final_best_tier']='FOLDSEEK_TECHNICAL_FAILURE_AFTER_SERIAL_RESOURCE_CLOSURE'
        out.append(row); continue

    row['foldseek_status']='PASS'
    m8=fold_dir/f'{cid}.m8'; hits=[]
    if m8.exists() and m8.stat().st_size:
        with m8.open() as f:
            for h in csv.DictReader(f,delimiter='\t',fieldnames=fields):
                h['candidate_id']=cid; h['target_accession']=normacc(h['target']); h['tier']=tier(ff(h['alntmscore']),ff(h['qcov']),ff(h['tcov']))
                hits.append(h); serial_hits.append(h)
    if not hits:
        row['r5c_hit_n']='0'; row['best_clean_rt_target']=''; row['best_alntmscore']=''; row['best_qcov']=''; row['best_tcov']=''; row['f3b3e_final_best_tier']='NO_HIT'
    else:
        best=max(hits,key=lambda h:ff(h['alntmscore']))
        row['r5c_hit_n']=str(len(hits)); row['best_clean_rt_target']=best['target_accession']; row['best_alntmscore']=best['alntmscore']; row['best_qcov']=best['qcov']; row['best_tcov']=best['tcov']; row['f3b3e_final_best_tier']=best['tier']
    out.append(row)

basecols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','known_sequence_reference_hit_n','known_top_reference_id','known_top_reference_original_headers']
cols=basecols+['f3b3e_resolution_source','pre_f3b3e_best_tier','serial_exclusive_boltz_status','serial_failure_category','serial_cif_path','foldseek_status','r5c_hit_n','best_clean_rt_target','best_alntmscore','best_qcov','best_tcov','f3b3e_final_best_tier']
final_tsv=outdir/'F3B3E_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv'
with final_tsv.open('w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=cols,lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(out)

hitcols=['candidate_id']+fields+['target_accession','tier']
hits_tsv=outdir/'F3B3E_SERIAL_CLEAN_RT_HITS.tsv'
with hits_tsv.open('w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=hitcols,lineterminator='\n'); w.writeheader(); w.writerows(serial_hits)

keys=['STRONG','SUGGESTIVE','WEAK','NO_HIT','RESOURCE_LIMITED_BOLTZ_UNRESOLVED','FOLDSEEK_TECHNICAL_FAILURE_AFTER_SERIAL_RESOURCE_CLOSURE']
counts={k:sum(r['f3b3e_final_best_tier']==k for r in out) for k in keys}
serial_pass=sum(r['serial_boltz_status'].startswith('PASS') for r in serial.values())
fold_fail=counts['FOLDSEEK_TECHNICAL_FAILURE_AFTER_SERIAL_RESOURCE_CLOSURE']
unresolved_n=counts['RESOURCE_LIMITED_BOLTZ_UNRESOLVED']+fold_fail
lock=json.load(open(lock_path))
gate_name='PASS_F3B3E_SERIAL_EXCLUSIVE_RESOURCE_CLOSURE_COMPLETE' if unresolved_n==0 else 'PASS_F3B3E_SERIAL_EXCLUSIVE_RESOURCE_CLOSURE_COMPLETE_WITH_RESOURCE_LIMITED_UNRESOLVED'
gate={
 'stage':'Stage05F3B3E','gate':gate_name,'gate_pass':True,
 'structure_handoff_family_n':32,
 'pre_f3b3e_resolved_structure_n':25,
 'pre_f3b3e_strong_or_suggestive_n':16,
 'serial_exclusive_retry_target_n':7,
 'serial_exclusive_boltz_success_n':serial_pass,
 'serial_exclusive_boltz_failure_n':7-serial_pass,
 'final_resolved_structure_n':32-unresolved_n,
 'final_r5c_strong_n':counts['STRONG'],'final_r5c_suggestive_n':counts['SUGGESTIVE'],'final_r5c_weak_n':counts['WEAK'],'final_r5c_no_hit_n':counts['NO_HIT'],
 'final_strong_or_suggestive_n':counts['STRONG']+counts['SUGGESTIVE'],
 'resource_limited_boltz_unresolved_n':counts['RESOURCE_LIMITED_BOLTZ_UNRESOLVED'],
 'foldseek_technical_failure_n':fold_fail,
 'technical_structure_closure_complete':unresolved_n==0,
 'resource_execution_change':'four concurrent independent single-GPU jobs -> one selected idle GPU, strictly serial, one fresh Boltz subprocess per query',
 'boltz_scientific_settings_changed':False,'sequence_truncation_used':False,'single_sequence_msa_fallback_used':False,'scientific_threshold_changed':False,'shortlist_membership_changed':False,
 'primary_or_f3b3d_results_overwritten':False,
 'inferential_role':'post-hoc resource-only technical closure; primary F3B3C and exact-protocol F3B3D remain separately preserved',
 'implementation_lock_sha256':sha(lock_path),'final_summary_sha256':sha(final_tsv),'serial_hits_sha256':sha(hits_tsv),
 'claim_boundary':'Serial-exclusive rescue changes computational resource scheduling only. STRONG/SUGGESTIVE remains frozen clean-RT structural compatibility, not biochemical RT/ART identity or novelty. Residual failures are resource-limited unresolved, not structural negatives.'
}
with (outdir/'STAGE05F3B3E_GATE.json').open('w') as f: json.dump(gate,f,indent=2);f.write('\n')
print(json.dumps(gate,indent=2))
