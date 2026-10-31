#!/usr/bin/env python3
import csv, hashlib, json, re, sys
from pathlib import Path

boltz_status=Path(sys.argv[1])
fold_status=Path(sys.argv[2])
fold_dir=Path(sys.argv[3])
audit=Path(sys.argv[4])
outdir=Path(sys.argv[5])
lock=Path(sys.argv[6])
outdir.mkdir(parents=True, exist_ok=True)

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

with boltz_status.open(newline='') as f:
    ss=list(csv.DictReader(f,delimiter='\t'))
assert len(ss)==32, len(ss)
with audit.open(newline='') as f:
    aa={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
with fold_status.open(newline='') as f:
    fs={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
L=json.load(open(lock))
fields=['query','target','evalue','bits','qcov','tcov','alnlen','qlen','tlen','fident','alntmscore']
out=[]; allhits=[]
for s in ss:
    cid=s['candidate_id'].strip(); bst=s['boltz_status'].strip(); cif=s['cif_path'].strip().rstrip('\r')
    base={
        'shortlist_rank':s['shortlist_rank'],'family_rep_id':s['family_rep_id'],'candidate_id':cid,
        'accession':s['accession'],'aa_length':s['aa_length'],'boltz_status':bst,'cif_path':cif,
        'foldseek_status':'NOT_RUN_BOLTZ_FAILURE',
        'known_sequence_reference_hit_n':aa.get(cid,{}).get('sequence_reference_hit_n',''),
        'known_top_reference_id':aa.get(cid,{}).get('top_reference_id',''),
        'known_top_reference_original_headers':aa.get(cid,{}).get('top_reference_original_headers','')
    }
    if not bst.startswith('PASS'):
        base.update({'r5c_hit_n':0,'best_clean_rt_target':'','best_alntmscore':'','best_qcov':'','best_tcov':'','best_tier':'BOLTZ_TECHNICAL_FAILURE'})
        out.append(base); continue
    fr=fs.get(cid)
    if fr is None or fr['foldseek_status']!='PASS':
        base['foldseek_status']=fr['foldseek_status'] if fr else 'MISSING_STATUS'
        base.update({'r5c_hit_n':0,'best_clean_rt_target':'','best_alntmscore':'','best_qcov':'','best_tcov':'','best_tier':'FOLDSEEK_TECHNICAL_FAILURE'})
        out.append(base); continue
    base['foldseek_status']='PASS'
    m8=fold_dir/f'{cid}.m8'; rows=[]
    if m8.exists() and m8.stat().st_size:
        with m8.open() as f:
            for r in csv.DictReader(f,delimiter='\t',fieldnames=fields):
                r['candidate_id']=cid; r['target_accession']=normacc(r['target']); r['tier']=tier(ff(r['alntmscore']),ff(r['qcov']),ff(r['tcov']))
                rows.append(r); allhits.append(r)
    if not rows:
        base.update({'r5c_hit_n':0,'best_clean_rt_target':'','best_alntmscore':'','best_qcov':'','best_tcov':'','best_tier':'NO_HIT'})
    else:
        best=max(rows,key=lambda r:ff(r['alntmscore']))
        base.update({'r5c_hit_n':len(rows),'best_clean_rt_target':best['target_accession'],'best_alntmscore':best['alntmscore'],'best_qcov':best['qcov'],'best_tcov':best['tcov'],'best_tier':best['tier']})
    out.append(base)

cols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','boltz_status','cif_path','foldseek_status','known_sequence_reference_hit_n','known_top_reference_id','known_top_reference_original_headers','r5c_hit_n','best_clean_rt_target','best_alntmscore','best_qcov','best_tcov','best_tier']
summary_tsv=outdir/'F3B3C_STRUCTURE_R5C_SUMMARY_32.tsv'
with summary_tsv.open('w',newline='\n') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=cols,lineterminator='\n'); w.writeheader(); w.writerows(out)
hit_tsv=outdir/'F3B3C_ALL_CLEAN_RT_HITS.tsv'
with hit_tsv.open('w',newline='\n') as f:
    hc=['candidate_id']+fields+['target_accession','tier']; w=csv.DictWriter(f,delimiter='\t',fieldnames=hc,lineterminator='\n'); w.writeheader(); w.writerows(allhits)
keys=['STRONG','SUGGESTIVE','WEAK','NO_HIT','BOLTZ_TECHNICAL_FAILURE','FOLDSEEK_TECHNICAL_FAILURE']
counts={k:sum(r['best_tier']==k for r in out) for k in keys}
complete=counts['BOLTZ_TECHNICAL_FAILURE']==0 and counts['FOLDSEEK_TECHNICAL_FAILURE']==0
gate={
 'stage':'Stage05F3B3C-v1.0.3',
 'gate':'PASS_F3B3C_FULL3D_AND_FROZEN_R5C_COMPLETE' if complete else 'PASS_F3B3C_PRIMARY_ATTEMPT_COMPLETE_WITH_TECHNICAL_FAILURES',
 'gate_pass':True,
 'hotfix':'CRLF_FOLDSEEK_HANDOFF_ONLY',
 'structure_handoff_family_n':32,
 'boltz_primary_attempt_n':len(out),
 'boltz_success_n':len(out)-counts['BOLTZ_TECHNICAL_FAILURE'],
 'boltz_technical_failure_n':counts['BOLTZ_TECHNICAL_FAILURE'],
 'foldseek_technical_failure_n':counts['FOLDSEEK_TECHNICAL_FAILURE'],
 'r5c_strong_n':counts['STRONG'],'r5c_suggestive_n':counts['SUGGESTIVE'],'r5c_weak_n':counts['WEAK'],'r5c_no_hit_n':counts['NO_HIT'],
 'r5c_strong_or_suggestive_survivor_n':counts['STRONG']+counts['SUGGESTIVE'],
 'structure_closure_complete':complete,
 'boltz_results_reused_from_v102':True,
 'boltz_recomputed':False,
 'boltz_model':'boltz2','msa_mode':'server-assisted','seed':20261007,'recycling_steps':3,'sampling_steps':200,'diffusion_samples':1,'no_kernels':True,
 'foldseek_evalue':1.0,'foldseek_max_seqs':100,'foldseek_exact_tmscore':1,
 'strong_rule':'alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50','suggestive_rule':'alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30',
 'no_single_sequence_fallback':True,'no_length_truncation':True,'shortlist_membership_changed':False,'scientific_threshold_changed':False,
 'known_matches_removed':False,'novelty_decision_performed':False,
 'implementation_lock_sha256':sha(lock),'summary_tsv_sha256':sha(summary_tsv),'all_hits_tsv_sha256':sha(hit_tsv),
 'claim_boundary':'STRONG/SUGGESTIVE denotes compatibility with the frozen clean-RT structural panel under the pre-existing R5C rule; it does not establish biochemical RT activity or ART identity. WEAK/NO_HIT is not evidence of non-RT or novelty. Boltz/Foldseek technical failures remain unresolved.'
}
with (outdir/'STAGE05F3B3C_GATE.json').open('w') as f: json.dump(gate,f,indent=2);f.write('\n')
print(json.dumps(gate,indent=2))
