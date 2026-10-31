#!/usr/bin/env python3
import csv, hashlib, json, shutil, sys
from pathlib import Path

gatep=Path(sys.argv[1]); tsvp=Path(sys.argv[2]); faap=Path(sys.argv[3]); outfa=Path(sys.argv[4]); outaud=Path(sys.argv[5])

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def fasta(path):
    hid=None; seq=[]
    with path.open() as f:
        for line in f:
            s=line.strip()
            if not s: continue
            if s.startswith('>'):
                if hid is not None: yield hid,''.join(seq)
                hid=s[1:].split()[0]; seq=[]
            else: seq.append(s)
        if hid is not None: yield hid,''.join(seq)

g=json.load(open(gatep))
assert g.get('gate_pass') is True
assert g.get('gate')=='PASS_F3V2_STAGE04B2_HISTORICAL_F3B3A_EXACT_RULE_PROJECTED_TO_DELTA_ONLY_FAMILIES'
assert g.get('structure_priority_shortlist_family_n')==32
assert g.get('known_rt_annotation_used_for_ranking') is False
assert g.get('novelty_used_for_ranking') is False
assert g.get('structure_used_for_ranking') is False
assert g.get('truth_used_for_ranking') is False
assert g.get('retuning_authorized') is False
assert g.get('scientific_threshold_changed') is False
assert sha(tsvp)==g['shortlist_tsv_sha256'], (sha(tsvp),g['shortlist_tsv_sha256'])
assert sha(faap)==g['shortlist_fasta_sha256'], (sha(faap),g['shortlist_fasta_sha256'])

with tsvp.open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
rows=[r for r in rows if str(r.get('shortlisted','')).lower() in ('true','1','yes')]
rows.sort(key=lambda r:int(r['shortlist_rank']))
assert len(rows)==32
assert [int(r['shortlist_rank']) for r in rows]==list(range(1,33))

seqs=dict(fasta(faap))
assert len(seqs)==32, len(seqs)

cols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length',
      'sequence_reference_hit_n','top_reference_id','top_reference_original_headers']
out=[]
seen=set()
for r in rows:
    cid=r['priority_candidate_id']
    acc=r['priority_candidate_accession']
    aa=int(r['priority_candidate_aa_length'])
    fam=r['family_rep_id']
    if cid in seen: raise SystemExit(f'FAIL duplicate priority candidate {cid}')
    seen.add(cid)
    if cid not in seqs: raise SystemExit(f'FAIL candidate absent from Stage04B2 FASTA: {cid}')
    if len(seqs[cid])!=aa: raise SystemExit(f'FAIL aa length mismatch {cid}: fasta={len(seqs[cid])} tsv={aa}')
    out.append({
      'shortlist_rank':r['shortlist_rank'],'family_rep_id':fam,'candidate_id':cid,
      'accession':acc,'aa_length':aa,
      'sequence_reference_hit_n':'NOT_EVALUATED_IN_DELTA_GENERALIZATION',
      'top_reference_id':'','top_reference_original_headers':''
    })

outfa.parent.mkdir(parents=True,exist_ok=True); outaud.parent.mkdir(parents=True,exist_ok=True)
shutil.copyfile(faap,outfa)
with outaud.open('w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=cols); w.writeheader(); w.writerows(out)

summary={
 'stage':'Stage05F3V2-Stage04B3-delta32-structure-handoff',
 'gate':'PASS_F3V2_STAGE04B3_DELTA32_HANDOFF_EXACT_FROM_STAGE04B2',
 'gate_pass':True,
 'source_stage04b2_gate_sha256':sha(gatep),
 'source_shortlist_tsv_sha256':sha(tsvp),
 'source_shortlist_fasta_sha256':sha(faap),
 'handoff_fasta_sha256':sha(outfa),
 'handoff_audit_sha256':sha(outaud),
 'handoff_family_n':32,
 'known_system_sequence_audit_performed':False,
 'shortlist_membership_changed':False,
 'science_change':False,
 'claim_boundary':'Exact Stage04B2 frozen 32-family compute-priority handoff to the historical F3B3C structure/R5C protocol. No known-system sequence annotation, novelty, structure, or truth information is used to alter membership.'
}
print(json.dumps(summary,indent=2))
