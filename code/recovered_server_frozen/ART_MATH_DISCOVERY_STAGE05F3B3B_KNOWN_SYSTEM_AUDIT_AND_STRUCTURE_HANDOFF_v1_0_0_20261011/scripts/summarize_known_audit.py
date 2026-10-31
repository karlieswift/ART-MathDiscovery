#!/usr/bin/env python3
import csv, hashlib, json, sys
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''):h.update(b)
    return h.hexdigest()

short_tsv,short_fa,ref_summary,ref_map,hits,outdir=sys.argv[1:7]
out=Path(outdir); out.mkdir(parents=True,exist_ok=True)

short=[]
with open(short_tsv) as f:
    r=csv.DictReader(f,delimiter='\t')
    for row in r: short.append(row)
assert len(short)==32
qids=[x['priority_candidate_id'] for x in short]
assert len(set(qids))==32

refsum=json.load(open(ref_summary))
refmeta={}
if Path(ref_map).is_file():
    with open(ref_map) as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r: refmeta[row['reference_id']]=row

byq={q:[] for q in qids}
if Path(hits).is_file() and Path(hits).stat().st_size:
    with open(hits) as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            q=row['query']
            if q not in byq: continue
            for k in ['fident','qcov','tcov','evalue','bits']:
                row[k]=float(row[k])
            row['alnlen']=int(row['alnlen'])
            byq[q].append(row)

rows=[]
for x in short:
    q=x['priority_candidate_id']; hs=byq[q]
    hs=sorted(hs,key=lambda z:(z['evalue'],-z['bits'],-z['qcov'],-z['fident'],z['target']))
    top=hs[0] if hs else None
    exact_any=any(h['fident']>=0.999999 and h['qcov']>=0.99 and h['tcov']>=0.99 for h in hs)
    high_any=any(h['fident']>=0.90 and h['qcov']>=0.80 and h['tcov']>=0.80 for h in hs)
    rec={
      'shortlist_rank':x['shortlist_rank'],'family_rep_id':x['family_rep_id'],'candidate_id':q,
      'accession':x['priority_candidate_accession'],'aa_length':x['priority_candidate_aa_length'],
      'best_rvt_profile':x['priority_candidate_best_formal_profile'],
      'sequence_reference_hit_n':len(hs),'exact_full_length_reference_match':exact_any,
      'high_identity_near_full_reference_match':high_any,
      'top_reference_id':top['target'] if top else '',
      'top_fident':top['fident'] if top else '', 'top_qcov':top['qcov'] if top else '',
      'top_tcov':top['tcov'] if top else '', 'top_evalue':top['evalue'] if top else '',
      'top_bits':top['bits'] if top else '',
      'top_reference_source_paths':refmeta.get(top['target'],{}).get('source_paths','') if top else '',
      'top_reference_original_headers':refmeta.get(top['target'],{}).get('original_headers','') if top else '',
    }
    rows.append(rec)

cols=list(rows[0])
with open(out/'F3B3B_KNOWN_SYSTEM_AUDIT_32.tsv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

# Handoff FASTA is copied by runner, byte identity checked there.
gate={
 'stage':'Stage05F3B3B',
 'gate':'PASS_F3B3B_KNOWN_SYSTEM_AUDIT_AND_STRUCTURE_HANDOFF',
 'gate_pass':True,
 'structure_handoff_family_n':32,
 'known_profile_annotation_complete':True,
 'known_sequence_reference_audit_performed':bool(refsum['sequence_reference_panel_available']),
 'reference_source_file_n':refsum['reference_source_file_n'],
 'reference_unique_aa_n':refsum['reference_unique_aa_n'],
 'candidate_with_sequence_reference_hit_n':sum(bool(byq[q]) for q in qids),
 'candidate_without_sequence_reference_hit_n':sum(not bool(byq[q]) for q in qids),
 'exact_full_length_reference_match_candidate_n':sum(any(h['fident']>=0.999999 and h['qcov']>=0.99 and h['tcov']>=0.99 for h in byq[q]) for q in qids),
 'high_identity_near_full_reference_match_candidate_n':sum(any(h['fident']>=0.90 and h['qcov']>=0.80 and h['tcov']>=0.80 for h in byq[q]) for q in qids),
 'known_matches_removed_from_structure_handoff':False,
 'shortlist_membership_changed':False,
 'novelty_decision_performed':False,
 'sequence_reference_negative_interpretation_allowed':False,
 'shortlist_tsv_sha256':sha(short_tsv),
 'shortlist_fasta_sha256':sha(short_fa),
 'audit_tsv_sha256':sha(out/'F3B3B_KNOWN_SYSTEM_AUDIT_32.tsv'),
 'claim_boundary':'Known-system profile/sequence similarity is descriptive post-freeze annotation only. Known matches are retained as positive calibration. Absence of a local sequence-reference hit is not evidence of novelty or non-RT identity.'
}
with open(out/'STAGE05F3B3B_GATE.json','w') as f:json.dump(gate,f,indent=2);f.write('\n')
print(json.dumps(gate,indent=2))
