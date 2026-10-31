#!/usr/bin/env python3
import argparse, hashlib, json
from pathlib import Path
import pandas as pd


def read_fasta(path):
    seqs={}; name=None; chunks=[]
    for raw in Path(path).read_text().splitlines():
        if raw.startswith('>'):
            if name is not None: seqs[name]=''.join(chunks)
            name=raw[1:].split()[0]
            chunks=[]
        else:
            chunks.append(raw.strip())
    if name is not None: seqs[name]=''.join(chunks)
    return seqs

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def boolish(x):
    return str(x).strip().lower() in {'true','1','yes'}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--d1',required=True); ap.add_argument('--d4',required=True)
    ap.add_argument('--r1',required=True); ap.add_argument('--r2',required=True)
    ap.add_argument('--out',default='outputs')
    a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)

    cand=pd.read_csv(Path(a.d1)/'config/NOVELTY_FIRST_999_CANDIDATES.tsv',sep='\t')
    cand=cand[cand['novelty_architecture_stratum'].astype(str).str.startswith(('N1_','N2_'))].copy()
    assert len(cand)==602
    cand['stratum_short']=cand['novelty_architecture_stratum'].str.extract(r'^(N[12])')[0]

    orf=pd.read_csv(Path(a.d1)/'config/NOVELTY_FIRST_24020_ORF_REGISTRY.tsv',sep='\t')
    links=pd.read_csv(Path(a.d1)/'config/NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv',sep='\t')
    uni=pd.read_csv(Path(a.d4)/'outputs/STAGE05D4_N1N2_PROTEIN_UNIVERSE.tsv',sep='\t')
    clu=pd.read_csv(Path(a.r1)/'outputs/STAGE05D4R1_NR90_CLUSTERS.tsv',sep='\t')
    tri=pd.read_csv(Path(a.r2)/'outputs/STAGE05D4R2_REMOTE_HOMOLOGY_TRIAGE.tsv',sep='\t')

    # normalize identifiers
    tri['representative_aa_sha256']=tri['query'].astype(str).str.replace(r'^AA_','',regex=True)
    assert tri['representative_aa_sha256'].nunique()==9245

    # Map every NR90 member AA -> all ORF instances -> N1/N2 candidate links.
    m=clu.merge(orf[['record_key','aa_sha256','aa_length','sequence_id','start0','end0','strand']],
                left_on='member_aa_sha256',right_on='aa_sha256',how='left',validate='m:m')
    if m['record_key'].isna().any():
        raise SystemExit(f"Missing ORF mappings for {int(m.record_key.isna().sum())} cluster-member rows")
    m=m.merge(links[['record_key','candidate_key','distance_to_math_candidate_bp']],on='record_key',how='left',validate='m:m')
    m=m.merge(cand[['candidate_key','region_rank','novelty_architecture_stratum','stratum_short']],on='candidate_key',how='inner',validate='m:1')
    if m.empty: raise SystemExit('No N1/N2 ORF-to-candidate mappings recovered')

    # representative-level deterministic architecture aggregation
    rows=[]
    for rep,g in m.groupby('representative_aa_sha256',sort=True):
        ur=uni[uni.aa_sha256.astype(str)==str(rep)]
        if len(ur)!=1:
            raise SystemExit(f"Representative {rep}: expected exactly one D4 universe row, got {len(ur)}")
        u=ur.iloc[0]
        distances=pd.to_numeric(g.distance_to_math_candidate_bp,errors='coerce')
        ranks=pd.to_numeric(g.region_rank,errors='coerce')
        keys=sorted(set(g.candidate_key.astype(str)))
        strata=sorted(set(g.stratum_short.astype(str)))
        rows.append({
          'representative_aa_sha256':rep,
          'query':'AA_'+rep,
          'aa_length':int(u.aa_length),
          'acidic_fraction':float(u.acidic_fraction),
          'gly_pro_fraction':float(u.gly_pro_fraction),
          'YxDD_present':boolish(u.YxDD_present),
          'D_x0_3_D_present':boolish(u.D_x0_3_D_present),
          'nr90_member_n':int(g.member_aa_sha256.nunique()),
          'orf_instance_n':int(g.record_key.nunique()),
          'candidate_key_n':len(keys),
          'arch_candidate_keys':';'.join(keys),
          'strata':';'.join(strata),
          'has_N1':bool('N1' in strata),
          'has_N2':bool('N2' in strata),
          'min_orf_to_candidate_distance_bp':int(distances.min()),
          'median_orf_to_candidate_distance_bp':float(distances.median()),
          'min_candidate_region_rank':int(ranks.min())
        })
    arch=pd.DataFrame(rows)
    z=tri.merge(arch,on=['query','representative_aa_sha256'],how='left',validate='1:1')
    if z['aa_length'].isna().any():
        raise SystemExit(f"Architecture mapping missing for {int(z.aa_length.isna().sum())} NR90 representatives")

    def tier(r):
        if r['triage_label']!='NO_SWISSPROT_HIT_NOT_NEGATIVE': return 'CURATED_SWISSPROT_HIT_TRIAGE'
        L=int(r['aa_length']); d=int(r['min_orf_to_candidate_distance_bp'])
        if L<180: return 'TIER_S_SHORT_LT180'
        if L>700: return 'TIER_L_LONG_GT700'
        if d<=1000: return 'TIER_A_COMPACT_STRUCTURE_READY'
        if d<=2500: return 'TIER_B_EXTENDED_STRUCTURE_READY'
        return 'TIER_C_REMOTE_DISTANCE_GT2500'
    z['structure_priority_tier']=z.apply(tier,axis=1)
    # transparent ordering, not a biological score
    tier_order={'TIER_A_COMPACT_STRUCTURE_READY':0,'TIER_B_EXTENDED_STRUCTURE_READY':1,
                'TIER_C_REMOTE_DISTANCE_GT2500':2,'TIER_S_SHORT_LT180':3,
                'TIER_L_LONG_GT700':4,'CURATED_SWISSPROT_HIT_TRIAGE':5}
    z['_tier_order']=z.structure_priority_tier.map(tier_order)
    z['_n1_order']=(~z.has_N1.astype(bool)).astype(int)
    z=z.sort_values(['_tier_order','_n1_order','min_orf_to_candidate_distance_bp','min_candidate_region_rank','representative_aa_sha256']).drop(columns=['_tier_order','_n1_order'])
    z.insert(0,'structure_priority_rank',range(1,len(z)+1))
    z.to_csv(out/'STAGE05D4R3_NR90_STRUCTURE_PRIORITY.tsv',sep='\t',index=False)

    # candidate-level audit: preserve all 602
    explode=[]
    for _,r in z.iterrows():
        for k in str(r.arch_candidate_keys).split(';'):
            explode.append((k,r.representative_aa_sha256,r.structure_priority_tier,int(r.min_orf_to_candidate_distance_bp)))
    ex=pd.DataFrame(explode,columns=['candidate_key','representative_aa_sha256','tier','rep_min_distance_bp'])
    cr=[]
    for _,c in cand.sort_values(['sequence_id','region_rank','candidate_key']).iterrows():
        g=ex[ex.candidate_key==c.candidate_key]
        cr.append({
          'candidate_key':c.candidate_key,'sequence_id':c.sequence_id,'region_rank':int(c.region_rank),
          'stratum':c.stratum_short,'novelty_architecture_stratum':c.novelty_architecture_stratum,
          'n_nr90_representatives':int(g.representative_aa_sha256.nunique()),
          'n_tier_A':int((g.tier=='TIER_A_COMPACT_STRUCTURE_READY').sum()),
          'n_tier_B':int((g.tier=='TIER_B_EXTENDED_STRUCTURE_READY').sum()),
          'has_tier_A':bool((g.tier=='TIER_A_COMPACT_STRUCTURE_READY').any()),
          'has_tier_A_or_B':bool(g.tier.isin(['TIER_A_COMPACT_STRUCTURE_READY','TIER_B_EXTENDED_STRUCTURE_READY']).any())
        })
    ctab=pd.DataFrame(cr)
    assert len(ctab)==602
    ctab.to_csv(out/'STAGE05D4R3_602_CANDIDATE_STRUCTURE_COVERAGE.tsv',sep='\t',index=False)

    # FASTA exports use frozen nonformal protein sequences; only representative IDs are emitted.
    seqs=read_fasta(Path(a.d4)/'outputs/STAGE05D4_N1N2_NONFORMAL_BROAD_RT_PROTEINS.faa')
    def write_tier(fname, tiers):
        q=z[z.structure_priority_tier.isin(tiers)]
        with open(out/fname,'w') as fh:
            for rep in q.representative_aa_sha256:
                key='AA_'+rep
                seq=seqs.get(key) or seqs.get(rep)
                if not seq: raise SystemExit(f"FASTA sequence missing for {key}")
                fh.write(f'>{key}\n')
                for i in range(0,len(seq),80): fh.write(seq[i:i+80]+'\n')
        return len(q)
    nA=write_tier('STAGE05D4R3_TIER_A_COMPACT_STRUCTURE_READY.faa',['TIER_A_COMPACT_STRUCTURE_READY'])
    nAB=write_tier('STAGE05D4R3_TIER_AB_STRUCTURE_READY.faa',['TIER_A_COMPACT_STRUCTURE_READY','TIER_B_EXTENDED_STRUCTURE_READY'])
    nAll=write_tier('STAGE05D4R3_ALL_NOHIT_180_700.faa',['TIER_A_COMPACT_STRUCTURE_READY','TIER_B_EXTENDED_STRUCTURE_READY','TIER_C_REMOTE_DISTANCE_GT2500'])

    counts=z.structure_priority_tier.value_counts().to_dict()
    summary={
      'n_nr90_queries':int(len(z)),
      'n_no_swissprot_hit':int((z.triage_label=='NO_SWISSPROT_HIT_NOT_NEGATIVE').sum()),
      'tier_counts':{k:int(v) for k,v in counts.items()},
      'n_tier_A_fasta':int(nA),'n_tier_AB_fasta':int(nAB),'n_all_nohit_180_700_fasta':int(nAll),
      'n_602_candidates':int(len(ctab)),
      'n_candidates_with_tier_A':int(ctab.has_tier_A.sum()),
      'n_candidates_with_tier_A_or_B':int(ctab.has_tier_A_or_B.sum())
    }
    (out/'STAGE05D4R3_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    gate={
      'stage':'Stage05D4-R3','gate':'PASS_COMPACT_STRUCTURE_QUERY_TIERS_FROZEN','gate_pass':True,
      'n_N1N2_candidates':602,'n_nr90_queries':int(len(z)),'n_no_swissprot_hit':summary['n_no_swissprot_hit'],
      'n_tier_A':int(nA),'n_tier_A_or_B':int(nAB),
      'dna_score_recomputed':False,'dna_rank_recomputed':False,'candidate_filtering_performed':False,
      'no_hit_interpreted_as_negative':False,'tier_membership_used_as_identity_assignment':False,
      'claim_boundary':'Tier A/B are frozen compute-priority sets for orthogonal structure/remote-homology analysis only; all 602 DNA-first N1/N2 candidates remain retained.'
    }
    (out/'GATE_STAGE05D4R3.json').write_text(json.dumps(gate,indent=2)+'\n')
    manifest=[]
    for p in sorted(out.glob('*')):
        if p.is_file() and p.name!='STAGE05D4R3_OUTPUT_MANIFEST.tsv': manifest.append((str(p),p.stat().st_size,sha(p)))
    pd.DataFrame(manifest,columns=['file','bytes','sha256']).to_csv(out/'STAGE05D4R3_OUTPUT_MANIFEST.tsv',sep='\t',index=False)
    print(json.dumps(gate,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
