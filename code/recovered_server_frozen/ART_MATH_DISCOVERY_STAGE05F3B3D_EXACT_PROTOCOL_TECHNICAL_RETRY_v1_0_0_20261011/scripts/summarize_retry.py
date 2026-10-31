#!/usr/bin/env python3
import csv, hashlib, json, re, sys
from pathlib import Path

primary_summary=Path(sys.argv[1])
retry_status_path=Path(sys.argv[2])
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

with primary_summary.open(newline='') as f:
    prim=list(csv.DictReader(f,delimiter='\t'))
with retry_status_path.open(newline='') as f:
    retry={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
with fold_status_path.open(newline='') as f:
    foldstat={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
assert len(prim)==32, len(prim)
assert len(retry)==10, len(retry)
fields=['query','target','evalue','bits','qcov','tcov','alnlen','qlen','tlen','fident','alntmscore']

out=[]; retry_hits=[]
for p in prim:
    row=dict(p)
    cid=p['candidate_id']
    primary_tier=p['best_tier']
    row['primary_best_tier']=primary_tier
    row['resolution_source']='PRIMARY'
    row['retry_boltz_status']='NOT_APPLICABLE'
    row['retry_cif_path']=''
    row['final_best_tier']=primary_tier
    if primary_tier!='BOLTZ_TECHNICAL_FAILURE':
        out.append(row); continue

    r=retry[cid]
    rbst=r['retry_boltz_status'].strip()
    row['retry_boltz_status']=rbst
    row['retry_cif_path']=r['cif_path'].strip().rstrip('\r')
    row['resolution_source']='EXACT_PROTOCOL_TECHNICAL_RETRY'
    if not rbst.startswith('PASS'):
        row['foldseek_status']='NOT_RUN_RETRY_BOLTZ_FAILURE'
        row['r5c_hit_n']='0'; row['best_clean_rt_target']=''; row['best_alntmscore']=''; row['best_qcov']=''; row['best_tcov']=''
        row['final_best_tier']='BOLTZ_TECHNICAL_FAILURE_AFTER_EXACT_RETRY'
        out.append(row); continue

    fs=foldstat.get(cid)
    if fs is None or fs['foldseek_status']!='PASS':
        row['foldseek_status']=fs['foldseek_status'] if fs else 'MISSING_STATUS'
        row['r5c_hit_n']='0'; row['best_clean_rt_target']=''; row['best_alntmscore']=''; row['best_qcov']=''; row['best_tcov']=''
        row['final_best_tier']='FOLDSEEK_TECHNICAL_FAILURE_AFTER_EXACT_RETRY'
        out.append(row); continue

    row['foldseek_status']='PASS'
    m8=fold_dir/f'{cid}.m8'; hits=[]
    if m8.exists() and m8.stat().st_size:
        with m8.open() as f:
            for h in csv.DictReader(f,delimiter='\t',fieldnames=fields):
                h['candidate_id']=cid; h['target_accession']=normacc(h['target']); h['tier']=tier(ff(h['alntmscore']),ff(h['qcov']),ff(h['tcov']))
                hits.append(h); retry_hits.append(h)
    if not hits:
        row['r5c_hit_n']='0'; row['best_clean_rt_target']=''; row['best_alntmscore']=''; row['best_qcov']=''; row['best_tcov']=''; row['final_best_tier']='NO_HIT'
    else:
        best=max(hits,key=lambda h:ff(h['alntmscore']))
        row['r5c_hit_n']=str(len(hits)); row['best_clean_rt_target']=best['target_accession']; row['best_alntmscore']=best['alntmscore']; row['best_qcov']=best['qcov']; row['best_tcov']=best['tcov']; row['final_best_tier']=best['tier']
    out.append(row)

# Stable output order and selected columns.
basecols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','known_sequence_reference_hit_n','known_top_reference_id','known_top_reference_original_headers']
cols=basecols+['resolution_source','primary_best_tier','retry_boltz_status','retry_cif_path','foldseek_status','r5c_hit_n','best_clean_rt_target','best_alntmscore','best_qcov','best_tcov','final_best_tier']
final_tsv=outdir/'F3B3D_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv'
with final_tsv.open('w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=cols,lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(out)

hitcols=['candidate_id']+fields+['target_accession','tier']
retry_hits_tsv=outdir/'F3B3D_RETRY_CLEAN_RT_HITS.tsv'
with retry_hits_tsv.open('w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=hitcols,lineterminator='\n'); w.writeheader(); w.writerows(retry_hits)

keys=['STRONG','SUGGESTIVE','WEAK','NO_HIT','BOLTZ_TECHNICAL_FAILURE_AFTER_EXACT_RETRY','FOLDSEEK_TECHNICAL_FAILURE_AFTER_EXACT_RETRY']
counts={k:sum(r['final_best_tier']==k for r in out) for k in keys}
retry_pass=sum(r['retry_boltz_status'].startswith('PASS') for r in retry.values())
unresolved=counts['BOLTZ_TECHNICAL_FAILURE_AFTER_EXACT_RETRY']+counts['FOLDSEEK_TECHNICAL_FAILURE_AFTER_EXACT_RETRY']
lock=json.load(open(lock_path))

gate={
 'stage':'Stage05F3B3D',
 'gate':'PASS_F3B3D_EXACT_PROTOCOL_TECHNICAL_RETRY_COMPLETE',
 'gate_pass':True,
 'primary_structure_handoff_n':32,
 'primary_boltz_success_n':22,
 'primary_boltz_technical_failure_n':10,
 'exact_protocol_retry_target_n':10,
 'exact_protocol_retry_boltz_success_n':retry_pass,
 'exact_protocol_retry_boltz_failure_n':10-retry_pass,
 'final_r5c_strong_n':counts['STRONG'],
 'final_r5c_suggestive_n':counts['SUGGESTIVE'],
 'final_r5c_weak_n':counts['WEAK'],
 'final_r5c_no_hit_n':counts['NO_HIT'],
 'final_strong_or_suggestive_n':counts['STRONG']+counts['SUGGESTIVE'],
 'unresolved_boltz_after_exact_retry_n':counts['BOLTZ_TECHNICAL_FAILURE_AFTER_EXACT_RETRY'],
 'unresolved_foldseek_after_exact_retry_n':counts['FOLDSEEK_TECHNICAL_FAILURE_AFTER_EXACT_RETRY'],
 'technical_structure_closure_complete':unresolved==0,
 'retry_rule':'all and only primary Boltz technical failures retried once under identical Boltz scientific settings',
 'diagnostic_category_used_for_retry_selection':False,
 'sequence_truncation_used':False,
 'single_sequence_msa_fallback_used':False,
 'scientific_threshold_changed':False,
 'shortlist_membership_changed':False,
 'primary_results_overwritten':False,
 'retry_results_are_separate_technical_rescue':True,
 'implementation_lock_sha256':sha(lock_path),
 'final_summary_sha256':sha(final_tsv),
 'retry_hits_sha256':sha(retry_hits_tsv),
 'claim_boundary':'Primary F3B3C remains the prespecified 32-query attempt. F3B3D is a separate exact-protocol technical retry of all primary Boltz failures. STRONG/SUGGESTIVE remains frozen clean-RT structural compatibility only, not biochemical RT/ART identity or novelty.'
}
with (outdir/'STAGE05F3B3D_GATE.json').open('w') as f:
    json.dump(gate,f,indent=2);f.write('\n')
print(json.dumps(gate,indent=2))
