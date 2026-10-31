#!/usr/bin/env python3
import argparse, hashlib, json
from pathlib import Path
from collections import defaultdict, Counter
import pandas as pd

TRANS=str.maketrans('ACGTN','TGCAN')

def canon_nt(s):
    s=str(s).upper().strip()
    rc=s.translate(TRANS)[::-1]
    return min(s,rc)

def split_semis(x):
    if pd.isna(x) or str(x).strip()=='': return []
    return [z for z in str(x).split(';') if z]

def mode_with_n(vals):
    vals=list(vals)
    if not vals: return ('',0)
    # deterministic lexical tiebreak
    c=Counter(vals)
    best=max(c.values())
    v=sorted([k for k,n in c.items() if n==best])[0]
    return v,best

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--d1-root',required=True)
    ap.add_argument('--e1-root',required=True)
    ap.add_argument('--fam50-tsv',required=True)
    ap.add_argument('--outdir',required=True)
    a=ap.parse_args()
    d1=Path(a.d1_root); e1=Path(a.e1_root); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    can=pd.read_csv(d1/'config/NOVELTY_FIRST_999_CANDIDATES.tsv',sep='\t')
    links=pd.read_csv(d1/'config/NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv',sep='\t')
    orf=pd.read_csv(d1/'config/NOVELTY_FIRST_24020_ORF_REGISTRY.tsv',sep='\t')
    rep=pd.read_csv(d1/'config/NOVELTY_FIRST_NONCODING_REPEAT_GROUPS.tsv',sep='\t')
    fam=pd.read_csv(a.fam50_tsv,sep='\t')
    ft=pd.read_csv(e1/'outputs/STAGE05E1_PRIMARY_REPLICATED_FAMILIES.tsv',sep='\t')
    ed=pd.read_csv(e1/'outputs/STAGE05E1_PRIMARY_EDGES.tsv',sep='\t')
    e1sum=json.loads((e1/'outputs/GATE_STAGE05E1.json').read_text())

    if not e1sum.get('gate_pass') or e1sum.get('primary_replicated_family_n')!=66 or e1sum.get('primary_replicated_candidate_n')!=324:
        raise SystemExit('Stage05E1 primary freeze mismatch')
    if len(can)!=999: raise SystemExit(f'Expected 999 D1 candidates, got {len(can)}')
    if len(ft)!=66: raise SystemExit(f'Expected 66 E1 primary families, got {len(ft)}')

    allkeys=set(can.candidate_key.astype(str))
    seqmap=can.set_index('candidate_key').sequence_id.astype(str).to_dict()

    # Full D1-999 qualifying repeat-description universe. Canonical core-string equality
    # is representation equality, not necessarily an exact biological repeat for N2.
    q=rep[(rep.candidate_key.astype(str).isin(allkeys)) &
          (rep.high_complexity==1) &
          (rep.full_span_outside_predicted_CDS==1) &
          (rep.overlap_candidate_bp>0)].copy()
    q['candidate_key']=q.candidate_key.astype(str)
    q['core_canonical']=q.core_nt.map(canon_nt)
    qe=q[q.assay.astype(str)=='EXACT_CORE_PRIMARY_DESCRIPTION'].copy()

    def prevalence(df):
        if len(df)==0: return {},{},{}
        cc=df.groupby('core_canonical').candidate_key.apply(lambda x:set(map(str,x))).to_dict()
        gg=df.groupby('core_canonical').sequence_id.apply(lambda x:set(map(str,x))).to_dict()
        ll=df.groupby('core_canonical').core_length.first().astype(int).to_dict()
        return cc,gg,ll
    any_core_cands,any_core_genomes,any_core_len=prevalence(q)
    exact_core_cands,exact_core_genomes,exact_core_len=prevalence(qe)

    # FAM50 membership from frozen D4R1; fan it back out to all D1 candidate-linked ORFs.
    fam_by_aa=defaultdict(set)
    for r in fam.itertuples(index=False):
        fam_by_aa[str(r.member_aa_sha256)].add(str(r.cluster_id))
    aa_by_record=orf.set_index('record_key').aa_sha256.astype(str).to_dict()
    fam_cands=defaultdict(set); fam_genomes=defaultdict(set)
    compact_rows=[]
    for r in links.itertuples(index=False):
        ck=str(r.candidate_key)
        if ck not in allkeys or int(r.distance_to_math_candidate_bp)>250: continue
        aa=aa_by_record.get(str(r.record_key),'')
        for fc in fam_by_aa.get(aa,()):
            fam_cands[fc].add(ck); fam_genomes[fc].add(seqmap[ck])
            compact_rows.append({'candidate_key':ck,'sequence_id':str(r.sequence_id),'record_key':str(r.record_key),
                                 'orf_start0':int(r.start0),'orf_end0':int(r.end0),'orf_strand':str(r.predicted_strand),
                                 'distance_to_math_candidate_bp':int(r.distance_to_math_candidate_bp),'aa_sha256':aa,'FAM50_cluster':fc})
    compact=pd.DataFrame(compact_rows)

    # Global prevalence tables for audit.
    core_rows=[]
    all_cores=sorted(set(any_core_cands)|set(exact_core_cands))
    for c in all_cores:
        core_rows.append({'core_canonical':c,'core_length':any_core_len.get(c,exact_core_len.get(c,0)),
                          'qualifying_candidate_n':len(any_core_cands.get(c,set())),'qualifying_genome_n':len(any_core_genomes.get(c,set())),
                          'exact_primary_candidate_n':len(exact_core_cands.get(c,set())),'exact_primary_genome_n':len(exact_core_genomes.get(c,set()))})
    pd.DataFrame(core_rows).to_csv(out/'STAGE05E2_GLOBAL_REPEAT_CORE_PREVALENCE.tsv',sep='\t',index=False)
    fprev=[]
    for fc in sorted(fam_cands):
        fprev.append({'FAM50_cluster':fc,'compact250_candidate_n':len(fam_cands[fc]),'compact250_genome_n':len(fam_genomes[fc])})
    pd.DataFrame(fprev).to_csv(out/'STAGE05E2_GLOBAL_FAM50_PREVALENCE.tsv',sep='\t',index=False)

    family_rows=[]; edge_rows=[]; signature={}
    for fr in ft.itertuples(index=False):
        fid=str(fr.family_id); members=set(split_semis(fr.candidate_keys))
        internal=ed[(ed.candidate_a.astype(str).isin(members))&(ed.candidate_b.astype(str).isin(members))].copy()
        cores=set(split_semis(fr.shared_repeat_cores_canonical)); fams=set(split_semis(fr.shared_FAM50_clusters))
        sig_any={c for c in cores if any_core_cands.get(c,set()) and any_core_cands.get(c,set())<=members and len(any_core_genomes.get(c,set()))>=2}
        sig_exact={c for c in cores if exact_core_cands.get(c,set()) and exact_core_cands.get(c,set())<=members and len(exact_core_genomes.get(c,set()))>=2}
        sig_fam={fc for fc in fams if fam_cands.get(fc,set()) and fam_cands.get(fc,set())<=members and len(fam_genomes.get(fc,set()))>=2}
        signature[fid]=(members,sig_any,sig_exact,sig_fam)
        joint_any=joint_exact=repeat_only=protein_only=broad=0
        for er in internal.itertuples(index=False):
            ecs=set(split_semis(er.shared_repeat_cores_canonical)); efs=set(split_semis(er.shared_FAM50_clusters))
            ha=bool(ecs & sig_any); he=bool(ecs & sig_exact); hf=bool(efs & sig_fam)
            if ha and hf: joint_any+=1
            if he and hf: joint_exact+=1
            if not (ha and hf):
                if ha: repeat_only+=1
                elif hf: protein_only+=1
                else: broad+=1
            edge_rows.append({'family_id':fid,'candidate_a':str(er.candidate_a),'candidate_b':str(er.candidate_b),
                              'shared_FAM50_cluster_n':int(er.shared_FAM50_cluster_n),'shared_repeat_core_string_n':int(er.shared_exact_repeat_core_n),
                              'has_family_exclusive_any_core':ha,'has_family_exclusive_exact_primary_core':he,'has_family_exclusive_compact_FAM50':hf,
                              'joint_any_signature':bool(ha and hf),'joint_exact_primary_signature':bool(he and hf),
                              'family_exclusive_any_cores':';'.join(sorted(ecs&sig_any)),
                              'family_exclusive_exact_primary_cores':';'.join(sorted(ecs&sig_exact)),
                              'family_exclusive_compact_FAM50':';'.join(sorted(efs&sig_fam))})
        if joint_exact>0: st='E2A_EXACT_PRIMARY_JOINT_SIGNATURE'
        elif joint_any>0: st='E2B_APPROX_OR_MIXED_JOINT_SIGNATURE'
        elif sig_any or sig_fam: st='E2C_PARTIAL_FAMILY_SPECIFIC_SIGNATURE'
        else: st='E2D_BROADLY_SHARED_EVIDENCE'
        family_rows.append({
            'family_id':fid,'candidate_n':int(fr.candidate_n),'genome_n':int(fr.genome_n),'N1_candidate_n':int(fr.N1_candidate_n),'N2_candidate_n':int(fr.N2_candidate_n),
            'min_region_rank':int(fr.min_region_rank),'edge_n':int(len(internal)),
            'family_exclusive_any_core_n':len(sig_any),'family_exclusive_exact_primary_core_n':len(sig_exact),'family_exclusive_compact_FAM50_n':len(sig_fam),
            'joint_any_signature_edge_n':joint_any,'joint_exact_primary_signature_edge_n':joint_exact,
            'joint_any_signature_edge_fraction':joint_any/len(internal) if len(internal) else 0.0,
            'joint_exact_primary_signature_edge_fraction':joint_exact/len(internal) if len(internal) else 0.0,
            'repeat_only_specific_edge_n':repeat_only,'protein_only_specific_edge_n':protein_only,'broad_evidence_edge_n':broad,
            'family_exclusive_any_cores':';'.join(sorted(sig_any)),'family_exclusive_exact_primary_cores':';'.join(sorted(sig_exact)),
            'family_exclusive_compact_FAM50_clusters':';'.join(sorted(sig_fam)),
            'stage05e2_stratum':st,'candidate_keys':str(fr.candidate_keys),'sequence_ids':str(fr.sequence_ids)
        })
    edge_audit=pd.DataFrame(edge_rows)
    edge_audit.to_csv(out/'STAGE05E2_EDGE_SIGNATURE_AUDIT.tsv',sep='\t',index=False)

    # Topology witness pairs. E2A uses exact-primary signature cores. E2B uses any qualifying signature core.
    witness=[]
    famdf=pd.DataFrame(family_rows)
    for frow in family_rows:
        fid=frow['family_id']; st=frow['stage05e2_stratum']; members,sig_any,sig_exact,sig_fam=signature[fid]
        if st=='E2A_EXACT_PRIMARY_JOINT_SIGNATURE': active_cores=sig_exact; rq=qe
        elif st=='E2B_APPROX_OR_MIXED_JOINT_SIGNATURE': active_cores=sig_any; rq=q
        else: continue
        if not active_cores or not sig_fam: continue
        for ck in sorted(members):
            reps=rq[(rq.candidate_key==ck)&(rq.core_canonical.isin(active_cores))]
            os=compact[(compact.candidate_key==ck)&(compact.FAM50_cluster.isin(sig_fam))]
            best=None
            for rp in reps.itertuples(index=False):
                for oo in os.itertuples(index=False):
                    if int(rp.end0)<=int(oo.orf_start0): gap=int(oo.orf_start0)-int(rp.end0); genomic='LEFT'
                    elif int(oo.orf_end0)<=int(rp.start0): gap=int(rp.start0)-int(oo.orf_end0); genomic='RIGHT'
                    else: gap=0; genomic='OVERLAP'
                    if genomic=='OVERLAP': tx='OVERLAP'
                    elif oo.orf_strand=='+': tx='UPSTREAM_5P' if genomic=='LEFT' else 'DOWNSTREAM_3P'
                    else: tx='DOWNSTREAM_3P' if genomic=='LEFT' else 'UPSTREAM_5P'
                    key=(gap,str(rp.core_canonical),str(oo.FAM50_cluster),str(oo.record_key),int(rp.start0),int(rp.end0))
                    row={'family_id':fid,'stage05e2_stratum':st,'candidate_key':ck,'sequence_id':seqmap[ck],
                         'repeat_group_id':str(rp.repeat_group_id),'repeat_assay':str(rp.assay),'signature_repeat_core':str(rp.core_canonical),
                         'signature_FAM50_cluster':str(oo.FAM50_cluster),'record_key':str(oo.record_key),'aa_sha256':str(oo.aa_sha256),
                         'repeat_start0':int(rp.start0),'repeat_end0':int(rp.end0),'orf_start0':int(oo.orf_start0),'orf_end0':int(oo.orf_end0),
                         'orf_strand':str(oo.orf_strand),'repeat_orf_gap_bp':gap,'genomic_repeat_side':genomic,'transcriptional_repeat_side':tx}
                    if best is None or key<best[0]: best=(key,row)
            if best is not None: witness.append(best[1])
    w=pd.DataFrame(witness)
    if len(w): w.to_csv(out/'STAGE05E2_SIGNATURE_WITNESS_PAIRS.tsv',sep='\t',index=False)
    else: pd.DataFrame(columns=['family_id']).to_csv(out/'STAGE05E2_SIGNATURE_WITNESS_PAIRS.tsv',sep='\t',index=False)

    # Merge topology summary onto family table.
    topo=[]
    if len(w):
        for fid,g in w.groupby('family_id'):
            mode,n=mode_with_n(g.transcriptional_repeat_side)
            topo.append({'family_id':fid,'signature_pair_candidate_n':int(g.candidate_key.nunique()),'signature_pair_genome_n':int(g.sequence_id.nunique()),
                         'median_repeat_orf_gap_bp':float(g.repeat_orf_gap_bp.median()),'max_repeat_orf_gap_bp':int(g.repeat_orf_gap_bp.max()),
                         'modal_transcriptional_repeat_side':mode,'modal_side_candidate_n':int(n),'modal_side_fraction':float(n/g.candidate_key.nunique())})
    topo=pd.DataFrame(topo)
    if len(topo): famdf=famdf.merge(topo,on='family_id',how='left')
    for c,default in [('signature_pair_candidate_n',0),('signature_pair_genome_n',0),('median_repeat_orf_gap_bp',float('nan')),('max_repeat_orf_gap_bp',float('nan')),('modal_transcriptional_repeat_side',''),('modal_side_candidate_n',0),('modal_side_fraction',float('nan'))]:
        if c not in famdf: famdf[c]=default

    order={'E2A_EXACT_PRIMARY_JOINT_SIGNATURE':0,'E2B_APPROX_OR_MIXED_JOINT_SIGNATURE':1,'E2C_PARTIAL_FAMILY_SPECIFIC_SIGNATURE':2,'E2D_BROADLY_SHARED_EVIDENCE':3}
    famdf['_order']=famdf.stage05e2_stratum.map(order)
    famdf=famdf.sort_values(['_order','genome_n','N1_candidate_n','joint_exact_primary_signature_edge_fraction','joint_any_signature_edge_fraction','candidate_n','family_id'],ascending=[True,False,False,False,False,False,True]).drop(columns='_order')
    famdf.insert(0,'review_order',range(1,len(famdf)+1))
    famdf.to_csv(out/'STAGE05E2_FAMILY_SPECIFICITY_TOPOLOGY.tsv',sep='\t',index=False)
    review_cols=['review_order','family_id','stage05e2_stratum','candidate_n','genome_n','N1_candidate_n','N2_candidate_n','min_region_rank',
                 'family_exclusive_exact_primary_core_n','family_exclusive_any_core_n','family_exclusive_compact_FAM50_n',
                 'joint_exact_primary_signature_edge_n','joint_any_signature_edge_n','joint_exact_primary_signature_edge_fraction','joint_any_signature_edge_fraction',
                 'signature_pair_candidate_n','signature_pair_genome_n','median_repeat_orf_gap_bp','modal_transcriptional_repeat_side','modal_side_fraction','candidate_keys','sequence_ids']
    famdf[review_cols].to_csv(out/'STAGE05E2_FAMILY_REVIEW_QUEUE.tsv',sep='\t',index=False)

    scount=Counter(famdf.stage05e2_stratum)
    candidate_by_stratum={}
    for st,g in famdf.groupby('stage05e2_stratum'): candidate_by_stratum[st]=int(g.candidate_n.sum())
    summary={
      'stage':'Stage05E2','n_D1_candidates_global_audit':999,'n_E1_primary_families':66,'n_E1_primary_candidates':324,
      'family_stratum_counts':dict(scount),'family_stratum_candidate_counts':candidate_by_stratum,
      'n_E2A_exact_primary_joint_families':int(scount.get('E2A_EXACT_PRIMARY_JOINT_SIGNATURE',0)),
      'n_E2B_approx_or_mixed_joint_families':int(scount.get('E2B_APPROX_OR_MIXED_JOINT_SIGNATURE',0)),
      'n_E2C_partial_specific_families':int(scount.get('E2C_PARTIAL_FAMILY_SPECIFIC_SIGNATURE',0)),
      'n_E2D_broad_shared_families':int(scount.get('E2D_BROADLY_SHARED_EVIDENCE',0)),
      'n_signature_witness_pairs':int(len(w)),
      'family_membership_changed':False,'candidate_filtering_performed':False,'uses_structure_results':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,
      'gate':'PASS_GLOBAL_SPECIFICITY_AND_TOPOLOGY_AUDIT_FROZEN','gate_pass':True,
      'claim_boundary':'E2A-E2D are review-order strata for replicated DNA+protein architectures. Family-exclusive means exclusive only within the frozen Stage05D1 999-candidate audit universe and does not establish RT activity, ncRNA function, or ART identity.'
    }
    (out/'STAGE05E2_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'GATE_STAGE05E2.json').write_text(json.dumps(summary,indent=2)+'\n')
    files=sorted(p for p in out.iterdir() if p.is_file() and p.name!='STAGE05E2_SHA256.tsv')
    with open(out/'STAGE05E2_SHA256.tsv','w') as h:
        for p in files: h.write(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
