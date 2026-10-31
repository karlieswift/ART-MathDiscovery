#!/usr/bin/env python3
import argparse,csv,hashlib,json,math,os
from pathlib import Path

def shafile(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--map',required=True)
    ap.add_argument('--cluster',required=True)
    ap.add_argument('--out',required=True)
    ap.add_argument('--mmseqs',required=True)
    ap.add_argument('--mmseqs-version',required=True)
    a=ap.parse_args(); OUT=Path(a.out); OUT.mkdir(parents=True,exist_ok=True)

    meta={}
    with open(a.map) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            meta[r['global_id']]=r
    fam={}
    seen_member={}
    membership=[]
    with open(a.cluster) as f:
        for line in f:
            line=line.rstrip('\n')
            if not line: continue
            rep,member=line.split('\t')[:2]
            if member not in meta: raise RuntimeError('Cluster member not in exact map: '+member)
            if member in seen_member: raise RuntimeError('Cluster member assigned multiple times: '+member)
            seen_member[member]=rep
            fam.setdefault(rep,[]).append(member)
    if set(seen_member)!=set(meta):
        miss=sorted(set(meta)-set(seen_member))[:20]
        extra=sorted(set(seen_member)-set(meta))[:20]
        raise RuntimeError(f'Cluster coverage mismatch missing={miss} extra={extra}')

    summary=[]
    classes={'OLD_ONLY':0,'DELTA_ONLY':0,'MIXED_OLD_DELTA':0}
    for rep,members in sorted(fam.items()):
        old=delta=shared=0
        for m in members:
            c=meta[m]['origin_class']
            if c in ('OLD_ONLY','SHARED_EXACT'): old+=1
            if c in ('DELTA_ONLY','SHARED_EXACT'): delta+=1
            if c=='SHARED_EXACT': shared+=1
            membership.append({'family_rep':rep,'member':m,'origin_class':c,'sequence_sha256':meta[m]['sequence_sha256']})
        cls='MIXED_OLD_DELTA' if old and delta else ('OLD_ONLY' if old else 'DELTA_ONLY')
        classes[cls]+=1
        summary.append({'family_rep':rep,'family_class':cls,'n_global_exact_sequences':len(members),'n_old_sequences':old,'n_delta_sequences':delta,'n_shared_exact_sequences':shared})

    with open(OUT/'COMBINED_FAMILY_MEMBERSHIP.tsv','w',newline='') as f:
        fields=['family_rep','member','origin_class','sequence_sha256']; w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(membership)
    with open(OUT/'COMBINED_FAMILY_SUMMARY.tsv','w',newline='') as f:
        fields=['family_rep','family_class','n_global_exact_sequences','n_old_sequences','n_delta_sequences','n_shared_exact_sequences']; w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(summary)

    pre=json.load(open(OUT/'PRECLUSTER_EXACT_CENSUS.json'))
    OLD_GB_NOMINAL=6.1216
    DELTA_GB=54.377819157
    old_formal=10598; delta_formal=88487
    old_unique=10329; delta_unique=85503
    old_families_frozen=4912; delta_families_frozen=26396
    census={
      'stage':'Stage05F3V2-Stage04A-60G-saturation-census',
      'gate':'PASS_F3V2_STAGE04A_60G_COMBINED_FAMILY_CENSUS_COMPLETE',
      'gate_pass':True,
      'old_universe_gb_nominal':OLD_GB_NOMINAL,
      'delta_universe_gb':DELTA_GB,
      'combined_universe_gb_nominal':OLD_GB_NOMINAL+DELTA_GB,
      'old_formal_rvt_n':old_formal,
      'delta_formal_rvt_n':delta_formal,
      'old_exact_unique_n':old_unique,
      'delta_exact_unique_n':delta_unique,
      'cross_universe_exact_overlap_n':pre['cross_universe_exact_overlap_n'],
      'combined_global_exact_unique_n':pre['combined_global_exact_unique_n'],
      'old_frozen_family_n':old_families_frozen,
      'delta_frozen_family_n':delta_families_frozen,
      'combined_reclustered_family_n':len(fam),
      'combined_old_only_family_n':classes['OLD_ONLY'],
      'combined_delta_only_family_n':classes['DELTA_ONLY'],
      'combined_mixed_old_delta_family_n':classes['MIXED_OLD_DELTA'],
      'delta_represented_combined_family_n':classes['DELTA_ONLY']+classes['MIXED_OLD_DELTA'],
      'old_represented_combined_family_n':classes['OLD_ONLY']+classes['MIXED_OLD_DELTA'],
      'delta_new_family_yield_per_Gb':classes['DELTA_ONLY']/DELTA_GB,
      'old_frozen_family_yield_per_Gb_nominal':old_families_frozen/OLD_GB_NOMINAL,
      'delta_frozen_family_yield_per_Gb':delta_families_frozen/DELTA_GB,
      'family_yield_density_ratio_delta_vs_old_nominal':(delta_families_frozen/DELTA_GB)/(old_families_frozen/OLD_GB_NOMINAL),
      'old_formal_rvt_yield_per_Gb_nominal':old_formal/OLD_GB_NOMINAL,
      'delta_formal_rvt_yield_per_Gb':delta_formal/DELTA_GB,
      'old_exact_unique_yield_per_Gb_nominal':old_unique/OLD_GB_NOMINAL,
      'delta_exact_unique_yield_per_Gb':delta_unique/DELTA_GB,
      'family_rule':'30% identity; 70% coverage; cov-mode0; cluster-mode0; sensitivity7.5',
      'mmseqs_path':a.mmseqs,
      'mmseqs_version':a.mmseqs_version,
      'global_exact_map_sha256':shafile(a.map),
      'combined_cluster_tsv_sha256':shafile(a.cluster),
      'science_change':False,
      'old_family_ids_redefined':False,
      'delta_family_ids_redefined':False,
      'retuning_authorized':False,
      'priority_reselection_performed':False,
      'claim_boundary':'Combined re-clustering is an analysis-only saturation census. It does not replace the frozen 6.12-Gb or 54.378-Gb family definitions and does not alter the five previously prioritized systems.'
    }
    with open(OUT/'STAGE04A_60G_SATURATION_CENSUS.json','w') as f: json.dump(census,f,indent=2); f.write('\n')
    print(json.dumps(census,indent=2))

if __name__=='__main__': main()
