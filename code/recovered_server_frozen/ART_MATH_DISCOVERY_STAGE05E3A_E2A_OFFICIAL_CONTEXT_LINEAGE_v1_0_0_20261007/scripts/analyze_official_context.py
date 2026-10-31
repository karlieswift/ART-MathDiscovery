#!/usr/bin/env python3
import argparse, json, hashlib, math, re
from pathlib import Path
from collections import Counter, defaultdict
import pandas as pd
from Bio import SeqIO

STRICT_RT=re.compile(r'(reverse transcriptase|RNA-directed DNA polymerase)',re.I)
PATTERNS={
 'retron':re.compile(r'\bretron\b',re.I),
 'dgr':re.compile(r'diversity[- ]generating|\bDGR\b',re.I),
 'crispr':re.compile(r'CRISPR|CRISPR-associated|\bCas\d*[A-Za-z-]*\b',re.I),
 'groupII':re.compile(r'group II intron|intron maturase',re.I),
 'telomerase':re.compile(r'telomerase',re.I),
 'retroelement':re.compile(r'retrotranspos|retroelement',re.I),
}

def qual1(f,k,default=''):
    v=f.qualifiers.get(k,[])
    return str(v[0]) if v else default

def parts_of_feature(f):
    return [(int(p.start),int(p.end)) for p in f.location.parts]

def interval_dist(a0,a1,b0,b1,L=None):
    def lin(x0,x1,y0,y1):
        if x1<=y0:return y0-x1
        if y1<=x0:return x0-y1
        return 0
    ds=[lin(a0,a1,b0,b1)]
    if L:
        ds += [lin(a0,a1,b0-L,b1-L),lin(a0,a1,b0+L,b1+L)]
    return min(ds)

def feature_dist(a0,a1,f,L=None):
    return min(interval_dist(a0,a1,b0,b1,L) for b0,b1 in parts_of_feature(f))

def overlap_len(a0,a1,f):
    return sum(max(0,min(a1,b1)-max(a0,b0)) for b0,b1 in parts_of_feature(f))

def feature_len(f):
    return sum(b-a for a,b in parts_of_feature(f))

def bounds(f):
    ps=parts_of_feature(f); return min(x[0] for x in ps),max(x[1] for x in ps)

def mode_frac(vals):
    vals=[x for x in vals if x not in ('',None) and not (isinstance(x,float) and math.isnan(x))]
    if not vals:return ('',0,0.0)
    c=Counter(vals); n=max(c.values()); v=sorted(k for k,z in c.items() if z==n)[0]
    return v,n,n/len(vals)

def context_flags(text):
    return {k:int(bool(rx.search(text))) for k,rx in PATTERNS.items()}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--frozen',required=True)
    ap.add_argument('--gbdir',required=True)
    ap.add_argument('--outdir',required=True)
    a=ap.parse_args(); fr=Path(a.frozen); gb=Path(a.gbdir); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    fam=pd.read_csv(fr/'STAGE05E3A_E2A_FAMILIES.tsv',sep='\t')
    cand=pd.read_csv(fr/'STAGE05E3A_E2A_CANDIDATES.tsv',sep='\t')
    wit=pd.read_csv(fr/'STAGE05E3A_E2A_WITNESS_PAIRS.tsv',sep='\t')
    freeze=json.loads((fr/'STAGE05E3A_FREEZE.json').read_text())
    if len(fam)!=12 or len(cand)!=119 or cand.sequence_id.nunique()!=61 or len(wit)!=80: raise SystemExit('Frozen Stage05E3A counts mismatch')
    wby=wit.set_index('candidate_key',drop=False).to_dict('index')
    records={}; meta=[]
    for sid in sorted(cand.sequence_id.unique()):
        p=gb/f'{sid}.gbff'
        if not p.exists(): raise SystemExit(f'Missing GenBank {p}')
        rs=list(SeqIO.parse(str(p),'genbank'))
        if not rs: raise SystemExit(f'Cannot parse {p}')
        rec=rs[0]; records[sid]=rec
        src=next((f for f in rec.features if f.type=='source'),None)
        tax=rec.annotations.get('taxonomy',[]) or []
        meta.append({'sequence_id':sid,'record_id':rec.id,'length_bp':len(rec.seq),'topology':rec.annotations.get('topology',''),
                     'organism':rec.annotations.get('organism',''),'taxonomy':' ; '.join(map(str,tax)),'taxonomy_tail3':' ; '.join(map(str,tax[-3:])),
                     'host':qual1(src,'host') if src else '', 'isolation_source':qual1(src,'isolation_source') if src else '',
                     'country':qual1(src,'country') if src else '', 'collection_date':qual1(src,'collection_date') if src else ''})
    pd.DataFrame(meta).to_csv(out/'STAGE05E3A_SEQUENCE_METADATA.tsv',sep='\t',index=False)
    md={r['sequence_id']:r for r in meta}
    nbh=[]; csum=[]; faa=[]; prot_seen=set()
    for cr in cand.itertuples(index=False):
        ck=str(cr.candidate_key); sid=str(cr.sequence_id); rec=records[sid]; L=len(rec.seq); circular=str(rec.annotations.get('topology','')).lower()=='circular'
        wr=wby.get(ck)
        if wr:
            a0=min(int(wr['repeat_start0']),int(wr['orf_start0'])); a1=max(int(wr['repeat_end0']),int(wr['orf_end0'])); anchor_source='E2A_EXACT_PRIMARY_JOINT_WITNESS'
        else:
            a0=int(cr.start0); a1=int(cr.end0); anchor_source='MATH_CANDIDATE_REGION_NO_FAMILY_EXCLUSIVE_WITNESS'
        cds=[f for f in rec.features if f.type=='CDS']
        local=[]
        for f in cds:
            d=feature_dist(a0,a1,f,L if circular else None)
            if d>10000: continue
            b0,b1=bounds(f); strand='+' if f.location.strand==1 else ('-' if f.location.strand==-1 else '')
            prod=qual1(f,'product'); pid=qual1(f,'protein_id'); gene=qual1(f,'gene'); locus=qual1(f,'locus_tag'); note=qual1(f,'note')
            text=' '.join([prod,gene,note])
            fl=context_flags(text); strict=int(bool(STRICT_RT.search(text)))
            row={'family_id':cr.family_id,'candidate_key':ck,'sequence_id':sid,'anchor_source':anchor_source,'anchor_start0':a0,'anchor_end0':a1,
                 'cds_start0':b0,'cds_end0':b1,'cds_strand':strand,'distance_to_anchor_bp':d,'within_5kb':int(d<=5000),
                 'protein_id':pid,'locus_tag':locus,'gene':gene,'product':prod,'note':note,'strict_rt_keyword':strict,**{f'flag_{k}':v for k,v in fl.items()}}
            local.append((f,row)); nbh.append(row)
            aa=qual1(f,'translation')
            if aa:
                key=hashlib.sha256(aa.encode()).hexdigest()
                if key not in prot_seen:
                    prot_seen.add(key); faa.append((key,sid,pid,locus,prod,aa))
        strict_n=sum(x[1]['strict_rt_keyword'] for x in local)
        fcounts={k:sum(x[1][f'flag_{k}'] for x in local) for k in PATTERNS}
        hypo=sum(1 for _,r in local if re.search(r'hypothetical|uncharacterized',str(r['product']),re.I))
        # map witness Prodigal ORF to official CDS
        match_class='NOT_APPLICABLE'; mpid=mprod=mgene=mlocus=''; mov=0.0; mbdelta=''; mstrand=''
        if wr:
            po0,po1=int(wr['orf_start0']),int(wr['orf_end0']); ps=str(wr['orf_strand']); best=None
            for f,row in local:
                ov=overlap_len(po0,po1,f)
                if ov<=0: continue
                den=min(max(1,po1-po0),max(1,feature_len(f))); frac=ov/den
                b0,b1=int(row['cds_start0']),int(row['cds_end0']); delta=abs(po0-b0)+abs(po1-b1); same=(row['cds_strand']==ps)
                key=(int(same),frac,-delta)
                if best is None or key>best[0]: best=(key,row,frac,delta)
            if best:
                _,br,mov,delta=best; mpid=br['protein_id']; mprod=br['product']; mgene=br['gene']; mlocus=br['locus_tag']; mstrand=br['cds_strand']; mbdelta=delta
                if mstrand==ps and delta<=30: match_class='SAME_STRAND_NEAR_EXACT_BOUNDARIES_15BP_EACH_SUM30'
                elif mstrand==ps and mov>=0.8: match_class='SAME_STRAND_HIGH_OVERLAP_GE80PCT'
                elif mov>=0.8: match_class='OPPOSITE_OR_UNKNOWN_STRAND_HIGH_OVERLAP_GE80PCT'
                else: match_class='PARTIAL_OVERLAP'
            else: match_class='NO_OFFICIAL_CDS_OVERLAP'
        csum.append({'family_id':cr.family_id,'review_order':cr.review_order,'candidate_key':ck,'sequence_id':sid,'region_rank':cr.region_rank,
                     'novelty_architecture_stratum':cr.novelty_architecture_stratum,'has_exact_primary_joint_witness':int(bool(wr)),
                     'anchor_source':anchor_source,'anchor_start0':a0,'anchor_end0':a1,'official_cds_10kb_n':len(local),'official_cds_5kb_n':sum(x[1]['within_5kb'] for x in local),
                     'hypothetical_or_uncharacterized_10kb_n':hypo,'strict_rt_keyword_10kb_n':strict_n,**{f'{k}_keyword_10kb_n':v for k,v in fcounts.items()},
                     'witness_official_match_class':match_class,'witness_official_overlap_fraction':mov,'witness_official_boundary_delta_sum_bp':mbdelta,
                     'witness_official_protein_id':mpid,'witness_official_locus_tag':mlocus,'witness_official_gene':mgene,'witness_official_product':mprod,
                     'organism':md[sid]['organism'],'taxonomy_tail3':md[sid]['taxonomy_tail3'],'host':md[sid]['host']})
    ndf=pd.DataFrame(nbh); ndf.to_csv(out/'STAGE05E3A_NEIGHBORHOOD_OFFICIAL_CDS.tsv',sep='\t',index=False)
    cdf=pd.DataFrame(csum); cdf.to_csv(out/'STAGE05E3A_CANDIDATE_OFFICIAL_CONTEXT.tsv',sep='\t',index=False)
    with open(out/'STAGE05E3A_NEIGHBORHOOD_OFFICIAL_PROTEINS.faa','w') as h:
        for sha,sid,pid,locus,prod,aa in faa:
            h.write(f'>{sha}|{sid}|{pid or locus}|{prod.replace(chr(9)," ")}\n')
            for i in range(0,len(aa),80): h.write(aa[i:i+80]+'\n')
    frows=[]
    fidx=fam.set_index('family_id')
    for fid,g in cdf.groupby('family_id'):
        base=fidx.loc[fid]
        wg=g[g.has_exact_primary_joint_witness==1]
        mm=wg.witness_official_match_class.astype(str)
        matched=mm.isin(['SAME_STRAND_NEAR_EXACT_BOUNDARIES_15BP_EACH_SUM30','SAME_STRAND_HIGH_OVERLAP_GE80PCT'])
        org_mode,org_n,org_frac=mode_frac(g.organism.astype(str))
        tax_mode,tax_n,tax_frac=mode_frac(g.taxonomy_tail3.astype(str))
        rt_c=int((g.strict_rt_keyword_10kb_n>0).sum())
        known_any=((g[[f'{k}_keyword_10kb_n' for k in PATTERNS]].sum(axis=1)>0)|(g.strict_rt_keyword_10kb_n>0)).sum()
        frows.append({'review_order':int(base.review_order),'family_id':fid,'candidate_n':int(base.candidate_n),'genome_n':int(base.genome_n),
                      'N1_candidate_n':int(base.N1_candidate_n),'N2_candidate_n':int(base.N2_candidate_n),
                      'signature_witness_candidate_n':int(len(wg)),'signature_witness_genome_n':int(wg.sequence_id.nunique()),
                      'official_same_strand_high_overlap_or_near_exact_witness_n':int(matched.sum()),
                      'official_same_strand_high_overlap_or_near_exact_witness_fraction':float(matched.mean()) if len(wg) else float('nan'),
                      'strict_rt_keyword_candidate_n_10kb':rt_c,'any_known_retro_system_keyword_candidate_n_10kb':int(known_any),
                      'unique_organism_name_n':int(g.organism.nunique()),'modal_organism':org_mode,'modal_organism_fraction':org_frac,
                      'unique_taxonomy_tail3_n':int(g.taxonomy_tail3.nunique()),'modal_taxonomy_tail3':tax_mode,'modal_taxonomy_tail3_fraction':tax_frac,
                      'median_repeat_orf_gap_bp_E2':base.median_repeat_orf_gap_bp,'modal_transcriptional_repeat_side_E2':base.modal_transcriptional_repeat_side,
                      'modal_side_fraction_E2':base.modal_side_fraction,
                      'context_interpretation':'KNOWN_RT_KEYWORD_CONTEXT_PRESENT' if rt_c>0 else ('OTHER_KNOWN_RETRO_SYSTEM_KEYWORD_CONTEXT_PRESENT' if known_any>0 else 'NO_CURATED_KNOWN_RETRO_SYSTEM_KEYWORD_IN_10KB_NOT_NOVELTY_PROOF')})
    fdf=pd.DataFrame(frows).sort_values('review_order'); fdf.to_csv(out/'STAGE05E3A_FAMILY_BIOLOGICAL_CONTEXT.tsv',sep='\t',index=False)
    summary={'stage':'Stage05E3A','n_E2A_families':12,'n_E2A_candidates':119,'n_sequence_ids':61,'n_signature_witness_candidates':80,
             'n_genbank_records_parsed':len(records),'n_neighborhood_official_cds_rows':int(len(ndf)),
             'n_unique_neighborhood_official_proteins':int(len(faa)),
             'n_witness_same_strand_high_overlap_or_near_exact_official_match':int(sum(cdf.witness_official_match_class.isin(['SAME_STRAND_NEAR_EXACT_BOUNDARIES_15BP_EACH_SUM30','SAME_STRAND_HIGH_OVERLAP_GE80PCT']))),
             'n_candidates_with_strict_rt_keyword_within_10kb':int((cdf.strict_rt_keyword_10kb_n>0).sum()),
             'n_families_with_strict_rt_keyword_within_10kb':int((fdf.strict_rt_keyword_candidate_n_10kb>0).sum()),
             'candidate_filtering_performed':False,'family_membership_changed':False,'selection_uses_structure_results':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,
             'gate':'PASS_E2A_OFFICIAL_CONTEXT_AND_LINEAGE_AUDIT_COMPLETE','gate_pass':len(records)==61,
             'claim_boundary':'Official GenBank annotation and taxonomy are post-discovery context. Keyword absence is not novelty evidence; keyword presence is triage evidence only. Predicted-to-official CDS overlap does not establish function.'}
    (out/'STAGE05E3A_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'GATE_STAGE05E3A.json').write_text(json.dumps(summary,indent=2)+'\n')
    files=sorted(p for p in out.iterdir() if p.is_file() and p.name!='STAGE05E3A_SHA256.tsv')
    with open(out/'STAGE05E3A_SHA256.tsv','w') as h:
        for p in files:h.write(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
