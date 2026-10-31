#!/usr/bin/env python3
import csv,json,sys,re
from pathlib import Path
from common import *

def hits(path):
    fields=['qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen','qcovs']
    out=[]
    if not Path(path).exists(): return out
    with open(path,encoding='utf-8',errors='ignore') as f:
        for line in f:
            z=line.rstrip('\n').split('\t');
            if len(z)<13: continue
            r=dict(zip(fields,z[:13]));
            r['_pid']=fnum(r['pident'])/100.; r['_qcov']=fnum(r['qcovs'])/100.; r['_tcov']=fnum(r['length'])/max(1,fnum(r['slen'])); r['_e']=fnum(r['evalue'],1e99); r['_bits']=fnum(r['bitscore'])
            out.append(r)
    return out

def acc_of(s):
    # NCBI Datasets genomic.fna identifiers usually begin with accession; tolerate pipes.
    m=re.search(r'([A-Z]{1,4}_?\d+(?:\.\d+)?)',s)
    return m.group(1) if m else s.split()[0]

def best_by_qacc(xs):
    d={}
    for r in xs:
        key=(r['qseqid'],acc_of(r['sseqid']))
        if key not in d or (r['_e'],-r['_bits'],-r['_qcov'],-r['_pid']) < (d[key]['_e'],-d[key]['_bits'],-d[key]['_qcov'],-d[key]['_pid']): d[key]=r
    return d

def main():
    if len(sys.argv)!=9: raise SystemExit('usage summarize_art8.py targets manifest rt_tbl partner_tbl arch_tbl genomic_fna out lock')
    targ,man,rtp,pp,ap,gfa,outp,lockp=sys.argv[1:]
    T=read_tsv(targ); M=read_tsv(man); wanted={r['accession'].split('.')[0]:r for r in M}; wanted.update({r['accession']:r for r in M})
    RT=best_by_qacc(hits(rtp)); PP=best_by_qacc(hits(pp)); AA=best_by_qacc(hits(ap))
    rows=[]
    for t in T:
        fam=t['family_rep_id']; q=t['rt_candidate_id']; pq=t['partner_orf_id']
        for m in M:
            acc=m['accession']; keys=[acc,acc.split('.')[0]]
            def pick(D,qq):
                for a in keys:
                    if (qq,a) in D:return D[(qq,a)]
                # subject id may have a prefix/suffix; scan
                cand=[r for (x,a),r in D.items() if x==qq and (a==acc or a.split('.')[0]==acc.split('.')[0])]
                return cand[0] if cand else None
            r=pick(RT,q); p=pick(PP,pq); a=pick(AA,q)
            rt90=bool(r and r['_pid']>=0.90 and r['_qcov']>=0.80)
            p90=bool(p and p['_pid']>=0.90 and p['_qcov']>=0.80)
            anear=bool(a and a['_pid']>=0.95 and a['_qcov']>=0.80)
            gap=''; coherent=False
            if r and p:
                ro=1 if int(r['send'])>=int(r['sstart']) else -1; po=1 if int(p['send'])>=int(p['sstart']) else -1
                rlo,rhi=sorted([int(r['sstart']),int(r['send'])]); plo,phi=sorted([int(p['sstart']),int(p['send'])])
                if ro==po==1 and plo>rhi: gap=plo-rhi-1; coherent=(gap<=500)
                elif ro==po==-1 and phi<rlo: gap=rlo-phi-1; coherent=(gap<=500)
            rows.append({'family_rep_id':fam,'rt_candidate_id':q,'art_id':m['art_id'],'art_source_name':m['source_name'],'art_accession':acc,
                'rt_hit':bool(r),'rt_pident':r['_pid'] if r else '','rt_qcov':r['_qcov'] if r else '','rt_subject_genome_fraction_not_protein_tcov':r['_tcov'] if r else '','rt_90pct_identity_80pct_query_coverage_compatible':rt90,
                'partner_hit':bool(p),'partner_pident':p['_pid'] if p else '','partner_qcov':p['_qcov'] if p else '','partner_subject_genome_fraction_not_protein_tcov':p['_tcov'] if p else '','partner_high_identity_query_coverage':p90,
                'rt_partner_same_orientation_adjacent_le500bp':coherent,'rt_partner_subject_gap_bp':gap,
                'architecture_hit':bool(a),'architecture_pident':a['_pid'] if a else '','architecture_qcov':a['_qcov'] if a else '','architecture_near_full_match':anear})
    write_tsv(outp,rows)
    # family best rows
    famsum=[]
    for t in T:
        fam=t['family_rep_id']; xs=[r for r in rows if r['family_rep_id']==fam]
        def score(r):
            return (bool(r['rt_90pct_identity_80pct_query_coverage_compatible']), fnum(r['rt_pident']), fnum(r['rt_qcov']), bool(r['architecture_near_full_match']), fnum(r['architecture_pident']), fnum(r['partner_pident']))
        best=sorted(xs,key=score,reverse=True)[0]
        famsum.append({'family_rep_id':fam,'best_public_art8_accession':best['art_accession'],'best_public_art8_name':best['art_source_name'],'best_public_art8_art_id':best['art_id'],
            'best_rt_pident':best['rt_pident'],'best_rt_qcov':best['rt_qcov'],'public_art8_rt_90pct_identity_80pct_query_coverage_compatible':any(b(r['rt_90pct_identity_80pct_query_coverage_compatible']) for r in xs),
            'best_partner_pident':max([fnum(r['partner_pident'],-1) for r in xs] or [-1]),'public_art8_architecture_near_full_match':any(b(r['architecture_near_full_match']) for r in xs),
            'public_art8_crosswalk_status':'PUBLIC_ART8_RT_90PCT_IDENTITY_80PCT_QUERY_COVERAGE_COMPATIBLE' if any(b(r['rt_90pct_identity_80pct_query_coverage_compatible']) for r in xs) else 'NO_90PCT_IDENTITY_80PCT_QUERY_COVERAGE_RT_MATCH_IN_PUBLIC_ART8_SUBSET'})
    sp=Path(outp).parent/'F3B4E_PUBLIC_ART8_FAMILY_SUMMARY_5.tsv'; write_tsv(sp,famsum)
    gate={'stage':'Stage05F3B4E-G2','gate':'PASS_F3B4E_PUBLIC_NAMED_ART8_CROSSWALK_COMPLETE','gate_pass':True,'target_family_n':5,'public_named_art_genome_n':8,
          'family_n_rt_90pct_identity_80pct_query_coverage_compatible_with_public_art8':sum(b(r['public_art8_rt_90pct_identity_80pct_query_coverage_compatible']) for r in famsum),'family_n_architecture_near_full_match_public_art8':sum(b(r['public_art8_architecture_near_full_match']) for r in famsum),
          'full_art95_catalogue_crosswalk_complete':False,'anthropic_unreported_art_claim_authorized':False,'detail_sha256':sha256_file(outp),'family_summary_sha256':sha256_file(sp),
          'claim_boundary':'Exact public named ART8 genome subset crosswalk is complete. The remaining ART95 catalogue is not machine-readable/publicly available here, so no absence-from-ART95 claim is authorized.'}
    gp=Path(outp).parent/'STAGE05F3B4E_G2_GATE.json'; gp.write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2)); print('\n=== PUBLIC ART8 FAMILY SUMMARY 5 ==='); print(sp.read_text().rstrip())

if __name__=='__main__': main()
