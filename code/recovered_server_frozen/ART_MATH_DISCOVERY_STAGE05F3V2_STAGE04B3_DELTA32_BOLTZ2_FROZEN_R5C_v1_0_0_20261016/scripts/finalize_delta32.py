#!/usr/bin/env python3
import csv, hashlib, json, shutil, sys
from pathlib import Path

hist_gate=Path(sys.argv[1]); hist_summary=Path(sys.argv[2]); hist_hits=Path(sys.argv[3]); stage04b2_gate=Path(sys.argv[4]); lockp=Path(sys.argv[5]); outdir=Path(sys.argv[6]); baseline_gate=Path(sys.argv[7])
outdir.mkdir(parents=True,exist_ok=True)

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

hg=json.load(open(hist_gate)); b2=json.load(open(stage04b2_gate)); L=json.load(open(lockp))
assert hg.get('gate_pass') is True
assert hg.get('structure_handoff_family_n')==32
assert b2.get('gate_pass') is True and b2.get('structure_priority_shortlist_family_n')==32
assert hg.get('shortlist_membership_changed') is False
assert hg.get('scientific_threshold_changed') is False

new_summary=outdir/'F3V2_STAGE04B3_STRUCTURE_R5C_SUMMARY_32.tsv'
new_hits=outdir/'F3V2_STAGE04B3_ALL_CLEAN_RT_HITS.tsv'
shutil.copyfile(hist_summary,new_summary); shutil.copyfile(hist_hits,new_hits)

with new_summary.open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==32
counts={k:sum(r['best_tier']==k for r in rows) for k in ['STRONG','SUGGESTIVE','WEAK','NO_HIT','TECHNICAL_FAILURE']}
completed=32-counts['TECHNICAL_FAILURE']
survivors=counts['STRONG']+counts['SUGGESTIVE']

gate={
 'stage':'Stage05F3V2-Stage04B3-delta32-Boltz2-frozen-R5C',
 'gate':'PASS_F3V2_STAGE04B3_DELTA32_FROZEN_STRUCTURE_R5C_COMPLETE' if counts['TECHNICAL_FAILURE']==0 else 'PASS_F3V2_STAGE04B3_DELTA32_PRIMARY_STRUCTURE_ATTEMPT_COMPLETE_WITH_TECHNICAL_FAILURES',
 'gate_pass':True,
 'stage04b2_gate_sha256':sha(stage04b2_gate),
 'implementation_lock_sha256':sha(lockp),
 'structure_handoff_family_n':32,
 'boltz_primary_attempt_n':32,
 'boltz_success_n':completed,
 'technical_failure_n':counts['TECHNICAL_FAILURE'],
 'r5c_strong_n':counts['STRONG'],
 'r5c_suggestive_n':counts['SUGGESTIVE'],
 'r5c_weak_n':counts['WEAK'],
 'r5c_no_hit_n':counts['NO_HIT'],
 'r5c_strong_or_suggestive_survivor_n':survivors,
 'survivor_fraction_of_all32':survivors/32.0,
 'survivor_fraction_of_resolved':(survivors/completed if completed else None),
 'structure_closure_complete':counts['TECHNICAL_FAILURE']==0,
 'boltz_model':'boltz2','msa_mode':'server-assisted','seed':20261007,'recycling_steps':3,'sampling_steps':200,'diffusion_samples':1,'no_kernels':True,
 'foldseek_evalue':1.0,'foldseek_max_seqs':100,'foldseek_exact_tmscore':1,
 'strong_rule':'alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50',
 'suggestive_rule':'alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30',
 'no_single_sequence_fallback':True,'no_length_truncation':True,
 'shortlist_membership_changed':False,'scientific_threshold_changed':False,
 'known_sequence_audit_used_for_handoff':False,'novelty_decision_performed':False,
 'historical_f3b3c_internal_gate_sha256':sha(hist_gate),
 'summary_tsv_sha256':sha(new_summary),'all_hits_tsv_sha256':sha(new_hits),
 'science_change':False,
 'claim_boundary':'This is a frozen out-of-universe structure generalization test of the Stage04B2 delta-only 32-family panel using the exact historical F3B3C Boltz-2 and clean-RT R5C settings. STRONG/SUGGESTIVE indicates structural compatibility with the clean-RT panel, not biochemical RT activity, ART identity, biological independence, or novelty. Technical failures remain unresolved.'
}
if baseline_gate.exists():
    try:
        bg=json.load(open(baseline_gate))
        if bg.get('gate_pass') is True:
            gate['historical_6p12gb_reference']={
              'gate_sha256':sha(baseline_gate),
              'boltz_success_n':bg.get('boltz_success_n'),
              'technical_failure_n':bg.get('technical_failure_n'),
              'r5c_strong_n':bg.get('r5c_strong_n'),
              'r5c_suggestive_n':bg.get('r5c_suggestive_n'),
              'r5c_strong_or_suggestive_survivor_n':bg.get('r5c_strong_or_suggestive_survivor_n'),
              'descriptive_only_not_used_for_delta_selection':True
            }
    except Exception:
        pass

outgate=outdir/'STAGE04B3_GATE.json'
with outgate.open('w') as f: json.dump(gate,f,indent=2); f.write('\n')
print(json.dumps(gate,indent=2))
