#!/usr/bin/env python3
import argparse, csv, hashlib, json
from pathlib import Path
from collections import defaultdict

OLD = {
    'bp': 6121620524,
    'formal_rvt_instance_n': 10598,
    'exact_aa_unique_n': 10329,
    'operational_family_n': 4912,
    'strict_plus_recurrent_partner_family_n': 7,
    'strong_clean_rt_family_n': 5,
}
DELTA_BP = 54377819157
HIST5 = {'F3BI37262050','F3BI37356863','F3BI27373723','F3BI55368798','F3BI55915933'}

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def read_json(p):
    with open(p) as f: return json.load(f)

def read_tsv(p):
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p, rows, cols):
    p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def discover(P, pattern, gate_rel, expected_gate):
    c=[]
    for d in sorted(P.glob(pattern)):
        p=d/gate_rel
        if not p.exists(): continue
        try: z=read_json(p)
        except Exception: continue
        if z.get('gate_pass') and z.get('gate')==expected_gate: c.append(d)
    if not c: raise SystemExit(f'HOLD no completed input root for {pattern}')
    return c[-1].resolve()

def splitset(s):
    return {x for x in (s or '').split(';') if x}

def closeout(a):
    P=Path(a.project_root).resolve(); H=Path(a.home).resolve(); O=H/'outputs'; O.mkdir(exist_ok=True)
    b5=discover(P,'ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B5_ALL_FORMAL_ARCHITECTURE_REPLAY_*','outputs/STAGE04B5_CLOSEOUT_GATE.json','PASS_STAGE04B5_54G_MATCHED_ARCHITECTURE_REPLAY_COMPLETE')
    b6=discover(P,'ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B6_ARCH41_CROSSWALK_AND_FRAGMENT_RESCUE_*','outputs/STAGE04B6_CLOSEOUT_GATE.json','PASS_STAGE04B6_ARCH41_CROSSWALK_AND_FRAGMENT_RESCUE_COMPLETE')
    b7=discover(P,'ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B7_RESCUED_DELTA_STRUCTURE_AND_REPLICATION_CLOSURE_*','outputs/STAGE04B7_CLOSEOUT_GATE.json','PASS_STAGE04B7_RESCUED_DELTA_STRUCTURE_AND_REPLICATION_CLOSURE')
    g5=read_json(b5/'outputs/STAGE04B5_CLOSEOUT_GATE.json')
    g6=read_json(b6/'outputs/STAGE04B6_CLOSEOUT_GATE.json')
    g7=read_json(b7/'outputs/STAGE04B7_CLOSEOUT_GATE.json')
    if (g5.get('formal_locus_n'),g5.get('exact_aa_unique_n'),g5.get('delta_operational_family_n'))!=(88487,85503,26396):
        raise SystemExit('HOLD B5 frozen census mismatch')
    arms={x['arm']:x for x in g5['architecture_by_arm']}
    mix=arms['MIXED_REPLICATION_ARM']; nov=arms['DELTA_ONLY_NOVELTY_ARM']
    if (mix['family_with_strict_plus_recurrent_partner_n'],nov['family_with_strict_plus_recurrent_partner_n'])!=(23,18):
        raise SystemExit('HOLD B5 architecture41 mismatch')
    rg=g6.get('rescue',{})
    if (rg.get('target_n'),rg.get('stop_bounded_both_n'),rg.get('edge_censored_unique_n'),rg.get('absent_n'))!=(18,12,1,5):
        raise SystemExit('HOLD B6 rescue mismatch')
    if (g7.get('primary_complete_structure_n'),g7.get('primary_positive_n'),g7.get('sensitivity_edge_censored_n'),g7.get('unresolved_absent_n'))!=(12,0,1,5):
        raise SystemExit('HOLD B7 final structure mismatch')
    # Semantic guard: connections are many-to-many, so do not call 5/5 system replication.
    cr=read_tsv(b7/'outputs/MIXED23_FAMILYWIDE_HISTORICAL_CROSSWALK.tsv')
    h5_to_mix=defaultdict(set); multi_mix=0; max_h5=0
    for r in cr:
        hs=splitset(r.get('historical5_overlap','')); max_h5=max(max_h5,len(hs))
        if len(hs)>1: multi_mix+=1
        for h in hs: h5_to_mix[h].add(r['delta_operational_family_id'])
    hist5_connected=set(h5_to_mix)&HIST5
    hist5_multi_mixed=sum(len(h5_to_mix[h])>1 for h in HIST5)
    if len(hist5_connected)!=5:
        raise SystemExit(f'HOLD expected 5 historical5 sequence-space connections, got {len(hist5_connected)}')
    many_to_many = multi_mix>0 or hist5_multi_mixed>0
    if not many_to_many: raise SystemExit('HOLD semantic guard expected many-to-many topology but none observed')

    arch41=mix['family_with_strict_plus_recurrent_partner_n']+nov['family_with_strict_plus_recurrent_partner_n']
    old_arch_frac=OLD['strict_plus_recurrent_partner_family_n']/OLD['operational_family_n']
    delta_arch_frac=arch41/g5['delta_operational_family_n']
    rows=[
      {'step':'frozen_selected_accessions','count':g5['selected_accession_n'],'role':'54G source universe'},
      {'step':'formal_rvt_profile_loci','count':g5['formal_locus_n'],'role':'profile-supported formal loci'},
      {'step':'operational_sequence_families','count':g5['delta_operational_family_n'],'role':'frozen delta family space'},
      {'step':'strict_plus_recurrent_partner_families_all','count':arch41,'role':'architecture-selected all delta families'},
      {'step':'mixed_architecture_families','count':mix['family_with_strict_plus_recurrent_partner_n'],'role':'old/delta mixed arm'},
      {'step':'delta_only_architecture_families','count':nov['family_with_strict_plus_recurrent_partner_n'],'role':'delta-only novelty arm'},
      {'step':'delta_only_complete_stop_bounded_primary','count':g7['primary_complete_structure_n'],'role':'primary structure-evaluable delta-only'},
      {'step':'delta_only_clean_rt_structural_positive','count':g7['primary_positive_n'],'role':'STRONG or SUGGESTIVE primary'},
      {'step':'delta_only_edge_censored_sensitivity','count':g7['sensitivity_edge_censored_n'],'role':'negative not authorized'},
      {'step':'delta_only_context_unresolved','count':g7['unresolved_absent_n'],'role':'technical/context unresolved'},
    ]
    write_tsv(O/'FINAL_54G_WATERFALL.tsv',rows,['step','count','role'])
    compare=[
      {'metric':'nucleotide_bp','historical_6p12G':OLD['bp'],'independent_54G':DELTA_BP,'interpretation':'searched DNA volume'},
      {'metric':'formal_rvt_profile_loci','historical_6p12G':OLD['formal_rvt_instance_n'],'independent_54G':g5['formal_locus_n'],'interpretation':'profile-level yield remains near proportional to volume'},
      {'metric':'operational_sequence_families','historical_6p12G':OLD['operational_family_n'],'independent_54G':g5['delta_operational_family_n'],'interpretation':'family-level novelty grows sublinearly'},
      {'metric':'strict_plus_recurrent_partner_families','historical_6p12G':OLD['strict_plus_recurrent_partner_family_n'],'independent_54G':arch41,'interpretation':'architecture endpoint recurs at similar family fraction'},
      {'metric':'strict_plus_recurrent_partner_family_fraction','historical_6p12G':f'{old_arch_frac:.8f}','independent_54G':f'{delta_arch_frac:.8f}','interpretation':f'delta/old ratio={delta_arch_frac/old_arch_frac:.4f}'},
      {'metric':'high_confidence_clean_rt_systems','historical_6p12G':OLD['strong_clean_rt_family_n'],'independent_54G':'0 of 12 evaluable delta-only architecture families','interpretation':'system-level novelty shows strong diminishing return under frozen closure'},
    ]
    write_tsv(O/'OLD6P12G_VS_DELTA54G_SUMMARY.tsv',compare,['metric','historical_6p12G','independent_54G','interpretation'])
    conn=[]
    for h in sorted(HIST5):
        fs=sorted(h5_to_mix[h])
        conn.append({'historical_family':h,'qualifying_mixed_delta_family_n':len(fs),'qualifying_mixed_delta_families':';'.join(fs),'authorized_label':'SEQUENCE_SPACE_CONNECTION_NOT_ONE_TO_ONE_SYSTEM_REPLICATION'})
    write_tsv(O/'HISTORICAL5_SEQUENCE_SPACE_CONNECTIONS.tsv',conn,list(conn[0].keys()))

    claim={
      'saturation':{
        'authorized': 'Strong diminishing-return / saturation signal for high-confidence system-level novelty within the frozen pipeline and sampled universes.',
        'not_authorized': [
          'The search space is exhausted.',
          'Another 100 Gb is unlikely to yield any new system.',
          'No additional ART-like architecture exists.'
        ],
        'reason': 'The 54G delta arm yielded 18 delta-only strict+recurrent architecture families but 0/12 complete stop-bounded primary ORFs passed frozen clean-RT structural thresholds; 5 families remain context unresolved and future universes may differ in ecology and assembly contiguity.'
      },
      'mixed_arm':{
        'historical5_with_sequence_space_connection_n':len(hist5_connected),
        'historical5_total_n':5,
        'mixed_family_with_multiple_historical5_connections_n':multi_mix,
        'historical5_connected_to_multiple_mixed_families_n':hist5_multi_mixed,
        'max_historical5_connections_for_one_mixed_family_n':max_h5,
        'sequence_space_representation_claim_authorized':True,
        'specific_one_to_one_system_replication_claim_authorized':False,
        'reason':'30% identity / 70% coverage operational connectivity is many-to-many.'
      },
      'anthropic_public_boundary':{
        'prioritized_6p12G_art_like_candidate_family_n':5,
        'public_art8_compatible_candidate_n':0,
        'full_art95_exclusion_complete':False,
        'absolute_global_sequence_or_locus_novelty_claim_authorized':False,
        'anthropic_unreported_art_claim_authorized':False,
        'authorized_narrative':'Independently identified ART-like architectures in public RT-containing loci that were not annotated as ART in the public records audited and had no compatible public ART8 counterpart.',
        'source':'Frozen prior Stage05F3B4E global public audit; not rerun in Stage04B8.'
      },
      'delta_only':{
        'strict_plus_recurrent_family_n':18,
        'primary_complete_stop_bounded_n':12,
        'primary_clean_rt_positive_n':0,
        'edge_censored_sensitivity_n':1,
        'context_unresolved_n':5,
        'new_computational_art_like_system_claim_authorized':False
      }
    }
    json.dump(claim,open(O/'FINAL_CLAIM_BOUNDARY.json','w'),indent=2)
    fig=[
      {'panel':'historical_6p12G','metric':'DNA_Gb','value':OLD['bp']/1e9},
      {'panel':'historical_6p12G','metric':'formal_RVT','value':OLD['formal_rvt_instance_n']},
      {'panel':'historical_6p12G','metric':'families','value':OLD['operational_family_n']},
      {'panel':'historical_6p12G','metric':'strict_plus_recurrent','value':7},
      {'panel':'historical_6p12G','metric':'STRONG_clean_RT_systems','value':5},
      {'panel':'independent_54G','metric':'DNA_Gb','value':DELTA_BP/1e9},
      {'panel':'independent_54G','metric':'formal_RVT','value':g5['formal_locus_n']},
      {'panel':'independent_54G','metric':'families','value':g5['delta_operational_family_n']},
      {'panel':'independent_54G','metric':'strict_plus_recurrent_all','value':arch41},
      {'panel':'independent_54G','metric':'strict_plus_recurrent_DELTA_ONLY','value':18},
      {'panel':'independent_54G','metric':'DELTA_ONLY_complete_primary','value':12},
      {'panel':'independent_54G','metric':'DELTA_ONLY_clean_RT_positive','value':0},
    ]
    write_tsv(O/'FIGURE6_FINAL_NUMBERS.tsv',fig,['panel','metric','value'])
    md=f'''# Stage04B8 manuscript-ready closeout\n\n## Final 54G result\nThe independent 54.378-Gb delta universe yielded 88,487 formal RVT-profile loci and 26,396 operational sequence families. Exhaustive architecture-first replay identified 41 families with strict upstream repeat geometry plus a recurrent downstream partner, comprising 23 MIXED and 18 DELTA_ONLY families. Among the 18 DELTA_ONLY families, 12 had uniquely recovered, complete stop-bounded ORFs suitable for primary structural interpretation; none passed the frozen clean-RT structural thresholds (0 STRONG, 0 SUGGESTIVE, 1 WEAK, 11 NO_HIT). One additional family was edge-censored and five remained context unresolved.\n\n## Saturation language\nAuthorized: **high-confidence system-level novelty shows a strong diminishing-return / saturation signal within the frozen pipeline and sampled universes.**\n\nNot authorized: **the search space is exhausted** or **another 100 Gb is unlikely to yield any new system**. A future orthogonal, contiguity-rich universe can still reveal rare systems.\n\n## Shared 6.12G / 54G pattern\nThe strict-repeat-plus-recurrent-partner family fraction is similar across the two universes: {OLD['strict_plus_recurrent_partner_family_n']}/{OLD['operational_family_n']} = {100*old_arch_frac:.4f}% historically versus {arch41}/{g5['delta_operational_family_n']} = {100*delta_arch_frac:.4f}% in the 54G delta universe. The key divergence occurs at orthogonal RT-fold closure: the historical architecture-selected set produced five STRONG clean-RT systems, whereas the DELTA_ONLY 54G primary set produced 0/12 structural positives.\n\n## Anthropic boundary\nThe five prioritized 6.12G candidates had no compatible public ART8 counterpart in the frozen public audit, but full ART95 exclusion was not complete. Therefore **Anthropic-unreported ART** and **absolute novel sequence/locus** claims remain unauthorized. The supported narrative is that the framework independently identified ART-like architectures in public RT-containing loci that were not annotated as ART in the public records audited.\n\n## MIXED-arm semantics\nAll five historical high-priority families have qualifying sequence-space connections to architecture-positive MIXED delta families. These connections are many-to-many under the operational 30% identity / 70% coverage rule, so they are evidence of sequence-space representation, not one-to-one independent system replication.\n'''
    (O/'MANUSCRIPT_READY_SUMMARY.md').write_text(md)
    prov=[]
    for role,p in [
      ('B5_gate',b5/'outputs/STAGE04B5_CLOSEOUT_GATE.json'),
      ('B6_gate',b6/'outputs/STAGE04B6_CLOSEOUT_GATE.json'),
      ('B7_gate',b7/'outputs/STAGE04B7_CLOSEOUT_GATE.json'),
      ('B7_mixed_crosswalk',b7/'outputs/MIXED23_FAMILYWIDE_HISTORICAL_CROSSWALK.tsv'),
      ('B7_structure_table',b7/'outputs/DELTA13_CORRECTED_STRUCTURAL_COMPATIBILITY.tsv')]:
        prov.append({'role':role,'path':str(p),'sha256':sha256(p)})
    write_tsv(O/'PROVENANCE_INPUTS.tsv',prov,['role','path','sha256'])
    gate={
      'stage':'Stage05F3V2-Stage04B8',
      'gate':'PASS_STAGE04B8_FINAL_54G_CLOSEOUT',
      'gate_pass':True,
      'formal_rvt_locus_n':g5['formal_locus_n'],
      'delta_operational_family_n':g5['delta_operational_family_n'],
      'strict_plus_recurrent_all_family_n':arch41,
      'mixed_strict_plus_recurrent_family_n':23,
      'delta_only_strict_plus_recurrent_family_n':18,
      'delta_only_primary_complete_n':12,
      'delta_only_primary_clean_rt_positive_n':0,
      'delta_only_edge_censored_n':1,
      'delta_only_context_unresolved_n':5,
      'historical5_with_sequence_space_connection_n':5,
      'specific_one_to_one_system_replication_claim_authorized':False,
      'anthropic_unreported_art_claim_authorized':False,
      'additional_100g_required_for_current_manuscript':False,
      'additional_100g_if_run':'OPTIONAL_INDEPENDENT_CONTIGUITY_AND_ECOLOGY_ORTHOGONAL_EXTENSION_ONLY',
      'new_discovery_run_performed':False,
      'structure_run_performed':False,
      'rules_changed':False,
      'retuning_performed':False,
      'final_54g_mainline_closed':True,
    }
    json.dump(gate,open(O/'STAGE04B8_FINAL_CLOSEOUT_GATE.json','w'),indent=2)
    with open(O/'SHA256SUMS.txt','w') as f:
        for p in sorted(x for x in O.iterdir() if x.is_file() and x.name!='SHA256SUMS.txt'):
            f.write(f'{sha256(p)}  {p.name}\n')
    print(json.dumps(gate,indent=2)); print('PASS_STAGE04B8_FINAL_54G_CLOSEOUT')

def selftest():
    old=7/4912; new=41/26396
    assert 0.0014 < old < 0.0015
    assert 0.0015 < new < 0.0016
    assert 0.9 < new/old < 1.2
    assert len(HIST5)==5
    print('PASS_STAGE04B8_SELFTEST')

def main():
    ap=argparse.ArgumentParser(); sp=ap.add_subparsers(dest='cmd',required=True)
    c=sp.add_parser('closeout'); c.add_argument('--project-root',required=True); c.add_argument('--home',required=True)
    sp.add_parser('selftest')
    a=ap.parse_args()
    if a.cmd=='closeout': closeout(a)
    else: selftest()
if __name__=='__main__': main()
