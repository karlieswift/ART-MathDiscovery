#!/usr/bin/env python3
import argparse, csv, hashlib, importlib.util, json, os, re, sys
from pathlib import Path

STOP={'TAA','TAG','TGA'}

def read_tsv(p):
    with open(p,newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p, rows, fields):
    p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,lineterminator='\n')
        w.writeheader()
        for r in rows:w.writerow(r)

def write_json(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w') as f:
        json.dump(x,f,indent=2,sort_keys=False); f.write('\n')

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def one_dir(P, token):
    ds=sorted([d for d in P.iterdir() if d.is_dir() and token in d.name])
    if not ds: raise SystemExit(f'HOLD cannot locate directory token={token}')
    return ds[-1]

def load_b6_module(B6):
    fp=B6/'scripts/stage04b6.py'
    if not fp.exists(): raise SystemExit(f'HOLD missing B6 code {fp}')
    spec=importlib.util.spec_from_file_location('b6mod',fp)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod

def original_contig_rescue(mod, genome, pc):
    # Orient the entire source contig so candidate translation is forward.
    if pc['strand']=='+':
        s=genome
        j=pc['start0']
        end=pc['end0']
    else:
        s=mod.rc(genome)
        j=len(genome)-pc['end0']
        end=len(genome)-pc['start0']
    cand=s[j:end]
    if len(cand)%3: raise RuntimeError('candidate nt not divisible by 3')
    # Scan in-frame to nearest stops exactly as the unitig rescue did.
    left=False; a=j; k=j-3
    while k>=0:
        if s[k:k+3] in STOP:
            a=k+3; left=True; break
        k-=3
    if not left: a=j%3
    right=False; b=end; k=end
    while k+3<=len(s):
        if s[k:k+3] in STOP:
            b=k; right=True; break
        k+=3
    if not right:
        b=len(s)-((len(s)-j)%3)
    aa=mod.tr(s[a:b])
    return {
        'original_contig_length_bp':len(s),
        'oriented_candidate_start0':j,
        'oriented_candidate_end0':end,
        'left_stop_bounded':left,
        'right_stop_bounded':right,
        'left_extension_aa':(j-a)//3,
        'right_extension_aa':(b-end)//3,
        'recovered_aa_length':len(aa),
        'recovered_aa':aa,
        'left_edge_distance_bp':j,
        'right_edge_distance_bp':len(s)-end,
    }

def discover_b3_counts(P):
    B3=one_dir(P,'STAGE04B3_DELTA32_BOLTZ2_FROZEN_R5C')
    # Prefer the frozen technical resource closure gate.
    gs=list(B3.glob('**/STAGE04B3_TECHNICAL_RESOURCE_CLOSURE_GATE.json'))
    for g in gs:
        try:x=json.load(open(g))
        except:continue
        c=x.get('closure',{})
        keys=['r5c_strong_n','r5c_suggestive_n','r5c_weak_n','r5c_no_hit_n']
        if all(k in c for k in keys):
            return B3,g,{
                'panel_n':32,
                'strong_n':int(c['r5c_strong_n']),
                'suggestive_n':int(c['r5c_suggestive_n']),
                'weak_n':int(c['r5c_weak_n']),
                'no_hit_n':int(c['r5c_no_hit_n']),
                'positive_n':int(c['r5c_strong_n'])+int(c['r5c_suggestive_n'])
            }
    # Fallback to summary TSV.
    ss=list(B3.glob('**/F3B3C_STRUCTURE_R5C_SUMMARY_32.tsv'))
    for s in ss:
        rows=read_tsv(s)
        if len(rows)!=32: continue
        tiercol='best_tier' if 'best_tier' in rows[0] else ('tier' if 'tier' in rows[0] else None)
        if not tiercol: continue
        cnt={k:sum(r[tiercol]==k for r in rows) for k in ['STRONG','SUGGESTIVE','WEAK','NO_HIT']}
        return B3,s,{'panel_n':32,'strong_n':cnt['STRONG'],'suggestive_n':cnt['SUGGESTIVE'],'weak_n':cnt['WEAK'],'no_hit_n':cnt['NO_HIT'],'positive_n':cnt['STRONG']+cnt['SUGGESTIVE']}
    raise SystemExit('HOLD cannot locate Stage04B3 32/32 structural closure counts')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('project_root')
    ap.add_argument('--outdir',default='outputs')
    a=ap.parse_args()
    P=Path(a.project_root).resolve(); O=Path(a.outdir).resolve(); O.mkdir(parents=True,exist_ok=True)

    B6=one_dir(P,'STAGE04B6_ARCH41_CROSSWALK_AND_FRAGMENT_RESCUE_v1_0_1')
    B7=one_dir(P,'STAGE04B7_RESCUED_DELTA_STRUCTURE_AND_REPLICATION_CLOSURE_v1_0_0')
    mod=load_b6_module(B6)

    b6tab=B6/'outputs/DELTA18_UNITIG_RESCUE.tsv'
    b7gate=B7/'outputs/STAGE04B7_CLOSEOUT_GATE.json'
    if not b6tab.exists() or not b7gate.exists(): raise SystemExit('HOLD B6/B7 outputs missing')
    b6=read_tsv(b6tab); b7=json.load(open(b7gate))
    absent=[r for r in b6 if r['rescue_status']=='ABSENT']
    if len(absent)!=5: raise SystemExit(f'HOLD expected 5 B6 ABSENT, got {len(absent)}')

    ledger=mod.discover_selected_ledger(P); lm=mod.ledger_map(ledger); source_root=ledger.parents[1]
    delta_faa=mod.discover_fasta(P,mod.EXPECTED_DELTA_FORMAL,('stage03','delta'))
    delta_aa=mod.fasta_dict(delta_faa)
    sf_cache={}; contig_cache={}; audit=[]
    for r in absent:
        cid=r['candidate_id']; pc=mod.parse_cid(cid); acc=pc['accession']
        if acc not in sf_cache:
            sf_cache[acc]=mod.find_source_file(source_root,acc,lm[acc],ledger)
        sf=sf_cache[acc]; key=(str(sf),pc['contig_id'])
        if key not in contig_cache:contig_cache[key]=mod.get_contig(sf,pc['contig_id'])
        genome=contig_cache[key]
        rr=original_contig_rescue(mod,genome,pc)
        aa0=delta_aa[cid].rstrip('*')
        # Candidate exact replay in original source.
        if pc['strand']=='+': cand_nt=genome[pc['start0']:pc['end0']]
        else: cand_nt=mod.rc(genome[pc['start0']:pc['end0']])
        aa=mod.tr(cand_nt).rstrip('*')
        if aa!=aa0: raise SystemExit(f'HOLD original-contig peptide replay {cid}')
        off=rr['left_extension_aa']
        if rr['recovered_aa'][off:off+len(aa)]!=aa:
            raise SystemExit(f'HOLD original-contig recovered peptide containment {cid}')
        complete=rr['left_stop_bounded'] and rr['right_stop_bounded']
        status='ORIGINAL_CONTIG_STOP_BOUNDED_BOTH' if complete else 'ORIGINAL_CONTIG_EDGE_CENSORED'
        audit.append({
            'delta_operational_family_id':r['delta_operational_family_id'],
            'candidate_id':cid,
            'accession':acc,
            'existing_aa_length':len(aa),
            'unitig_status':'ABSENT',
            'source_contig_file':str(sf),
            'source_contig_length_bp':rr['original_contig_length_bp'],
            'left_edge_distance_bp':rr['left_edge_distance_bp'],
            'right_edge_distance_bp':rr['right_edge_distance_bp'],
            'left_stop_bounded':rr['left_stop_bounded'],
            'right_stop_bounded':rr['right_stop_bounded'],
            'original_contig_recovered_aa_length':rr['recovered_aa_length'],
            'left_extension_aa':rr['left_extension_aa'],
            'right_extension_aa':rr['right_extension_aa'],
            'original_contig_status':status,
            'structure_followup_role':('POST_B8_TECHNICAL_COMPLETE_ORF_SENSITIVITY' if complete else 'KEEP_CONTEXT_UNRESOLVED')
        })
    fields=list(audit[0].keys())
    write_tsv(O/'UNRESOLVED5_ORIGINAL_CONTIG_STOP_BOUNDARY_AUDIT.tsv',audit,fields)

    B3,b3source,b3=discover_b3_counts(P)
    if b3['panel_n']!=32: raise SystemExit('HOLD B3 panel n')
    # Frozen project result expected from Stage04B3 technical closure.
    # We do not hard-code the counts; we read them from the server artifact.
    arch={
        'delta_only_strict_plus_recurrent_family_n':18,
        'primary_complete_stop_bounded_n':int(b7['primary_complete_structure_n']),
        'primary_strong_n':int(b7['primary_strong_n']),
        'primary_suggestive_n':int(b7['primary_suggestive_n']),
        'primary_positive_n':int(b7['primary_positive_n']),
        'edge_censored_n':int(b7['sensitivity_edge_censored_n']),
        'unitig_absent_n':5
    }

    complete5=sum(r['original_contig_status']=='ORIGINAL_CONTIG_STOP_BOUNDED_BOTH' for r in audit)
    reconcile={
        'stage':'Stage05F3V2-Stage04B8R1',
        'gate':'PASS_STAGE04B8R1_CROSS_ROUTE_RECONCILIATION_AND_UNRESOLVED5_AUDIT',
        'gate_pass':True,
        'protein_priority_route':{
            'source_space':'24,810 DELTA_ONLY combined families',
            'selection':'historical F3B3A frozen 32-family projection: replication/length/model-coverage/motif-topology',
            'stage04b3_source':str(b3source),
            **b3
        },
        'architecture_first_route':{
            'source_space':'all 88,487 formal delta loci -> strict upstream repeat + recurrent partner -> DELTA_ONLY',
            **arch
        },
        'unresolved5_original_contig_audit':{
            'unitig_absent_n':5,
            'original_contig_stop_bounded_both_n':complete5,
            'original_contig_edge_censored_n':5-complete5,
            'table':str((O/'UNRESOLVED5_ORIGINAL_CONTIG_STOP_BOUNDARY_AUDIT.tsv').resolve())
        },
        'claim_correction':{
            '54g_has_clean_rt_structural_positives': b3['positive_n']>0,
            '54g_overall_clean_rt_absence_claim_authorized':False,
            'high_confidence_system_level_saturation_claim_authorized':False,
            'reason':'The 54G protein-priority DELTA_ONLY panel contains clean-RT structural positives, whereas the separate architecture-first DELTA_ONLY set has 0 structural positives among its currently primary-evaluable ORFs. These are different frozen panels answering different questions. Protein-positive loci were previously architecture-context censored in fragmented assemblies, so lack of joint closure cannot be interpreted as exhaustion.',
            'authorized_language':'In the 54G delta universe, frozen protein-priority selection retained RT-fold generalization, while exhaustive architecture-first selection exposed additional strict+recurrent families whose evaluable short ORFs lacked clean-RT fold support. The two axes did not achieve joint system-level closure in the fragmented delta assemblies.'
        },
        'science_change':False,
        'retuning_performed':False,
        'new_discovery_performed':False,
        'structure_run_performed':False
    }
    write_json(O/'STAGE04B8R1_RECONCILIATION_GATE.json',reconcile)
    with open(O/'MANUSCRIPT_CORRECTION.md','w') as f:
        f.write('# Stage04B8R1 semantic correction\n\n')
        f.write('Stage04B8 v1.0.0 branch-specific counts remain valid, but its manuscript-ready saturation sentence must not be read as a statement that the 54G universe contained no clean-RT structural positives overall.\n\n')
        f.write(f'The protein-priority DELTA32 route yielded {b3["positive_n"]}/32 clean-RT structural positives ({b3["strong_n"]} STRONG, {b3["suggestive_n"]} SUGGESTIVE), whereas the separate architecture-first DELTA_ONLY route yielded {arch["primary_positive_n"]}/{arch["primary_complete_stop_bounded_n"]} positives among currently primary-complete ORFs.\n\n')
        f.write('Therefore the correct interpretation is cross-axis non-closure under severe assembly-context censoring, not proof of system-space saturation.\n')
    # provenance hashes
    with open(O/'SHA256SUMS.txt','w') as f:
        for p in sorted(x for x in O.iterdir() if x.is_file() and x.name!='SHA256SUMS.txt'):
            f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(reconcile,indent=2))
    print('PASS_STAGE04B8R1_CROSS_ROUTE_RECONCILIATION_AND_UNRESOLVED5_AUDIT')

if __name__=='__main__':main()
