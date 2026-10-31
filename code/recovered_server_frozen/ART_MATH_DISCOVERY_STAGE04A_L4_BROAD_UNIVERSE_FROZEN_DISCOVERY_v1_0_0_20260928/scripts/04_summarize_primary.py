#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, json, math, statistics
from common import sha256_file, write_json

def qtile(xs,p):
    if not xs: return None
    ys=sorted(xs); pos=(len(ys)-1)*p; lo=int(math.floor(pos)); hi=int(math.ceil(pos));
    return ys[lo] if lo==hi else ys[lo]+(ys[hi]-ys[lo])*(pos-lo)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); a=ap.parse_args(); root=Path(a.root).resolve()
    with open(root/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv') as f: cohort=list(csv.DictReader(f,delimiter='\t'))
    byid={r['analysis_id']:r for r in cohort}; per=[]; allregs=[]
    for aid,r in byid.items():
        gdir=root/'outputs/genomes'/aid; gate=json.loads((gdir/'STAGE03A_PRIMARY_GATE.json').read_text()); rp=gdir/'STAGE03A_UNION_REGIONS.tsv'
        with open(rp) as f: regs=list(csv.DictReader(f,delimiter='\t'))
        L=int(r['length_bp']); n=len(regs); burden=n/(L/1e6)
        fam_counts={'ACNV_FAMILY':0,'IRREGULARITY_V3_D20':0,'GENERIC_RECURRENCE':0}
        for x in regs:
            fams=set(filter(None,x.get('support_families','').split(';')))
            for fam in fam_counts: fam_counts[fam]+=int(fam in fams)
            y={'analysis_id':aid,'sequence_id':r['sequence_id'],'arm':r['arm'],'length_bp':L,**x}; allregs.append(y)
        rank1=regs[0] if regs else {}
        per.append({'analysis_id':aid,'sequence_id':r['sequence_id'],'arm':r['arm'],'source':r['source'],'release_date':r['release_date'],'cluster_id':r.get('cluster_id',''),'length_bp':L,'n_union_regions':n,'candidate_burden_per_Mb':burden,'has_any_union_region':int(n>0),'rank1_region_id':rank1.get('region_id',''),'rank1_start0':rank1.get('start0',''),'rank1_end0':rank1.get('end0',''),'rank1_support_family_n':rank1.get('support_family_n',''),'rank1_support_families':rank1.get('support_families',''),'n_regions_with_ACNV':fam_counts['ACNV_FAMILY'],'n_regions_with_IRR':fam_counts['IRREGULARITY_V3_D20'],'n_regions_with_RECURRENCE':fam_counts['GENERIC_RECURRENCE'],'stage03a_primary_gate':gate.get('gate','')})
    ppath=root/'outputs/STAGE04A_PER_SEQUENCE_BURDEN.tsv'; fields=list(per[0].keys()) if per else []
    with open(ppath,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(per)
    arpath=root/'outputs/STAGE04A_ALL_FROZEN_UNION_REGIONS.tsv'
    if allregs:
        fields2=list(allregs[0].keys());
        with open(arpath,'w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields2,delimiter='\t'); w.writeheader(); w.writerows(allregs)
    else: arpath.write_text('analysis_id\tsequence_id\tarm\n')
    arms={}
    for arm in sorted(set(r['arm'] for r in per)):
        z=[r for r in per if r['arm']==arm]; burdens=[float(r['candidate_burden_per_Mb']) for r in z]; totalbp=sum(int(r['length_bp']) for r in z); totalreg=sum(int(r['n_union_regions']) for r in z)
        d={'n_sequences':len(z),'total_Mb':totalbp/1e6,'total_union_regions':totalreg,'pooled_candidate_burden_per_Mb':totalreg/(totalbp/1e6) if totalbp else None,'median_per_sequence_burden_per_Mb':statistics.median(burdens) if burdens else None,'iqr_per_sequence_burden_per_Mb':[qtile(burdens,.25),qtile(burdens,.75)] if burdens else None,'p95_per_sequence_burden_per_Mb':qtile(burdens,.95),'fraction_sequences_ge1':sum(int(r['has_any_union_region']) for r in z)/len(z) if z else None}
        clusters={}
        for r in z:
            cid=(r.get('cluster_id') or '').strip()
            if cid: clusters[cid]=clusters.get(cid,0)+1
        if clusters and all((r.get('cluster_id') or '').strip() for r in z):
            wreg=wbp=0.0
            for r in z:
                w=1.0/clusters[r['cluster_id'].strip()]; wreg+=w*int(r['n_union_regions']); wbp+=w*int(r['length_bp'])
            d['cluster_weighted_candidate_burden_per_Mb']=wreg/(wbp/1e6) if wbp else None; d['n_clusters']=len(clusters)
        else:
            d['cluster_weighted_candidate_burden_per_Mb']=None; d['n_clusters']=None
        arms[arm]=d
    summary={'stage':'Stage04A','gate':'PASS_STAGE04A_FROZEN_ENGINE_BROAD_UNIVERSE_PRIMARY_SCORES_FROZEN','gate_pass':True,'truth_used':False,'annotation_used_for_scoring':False,'top_k_survival_cutoff_used':False,'retuning_authorized':False,'n_sequences':len(per),'n_total_union_regions':len(allregs),'arms':arms,'interpretation':'Candidate burden and frozen priority ranks only. Rank 1 is a mathematical anomaly priority, not an ART identity call. Novelty/RT/partner annotation is deferred to Stage04B.'}
    write_json(root/'outputs/STAGE04A_COHORT_SUMMARY.json',summary)
    gate={**summary,'universe_manifest_sha256':sha256_file(root/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv'),'per_sequence_burden_sha256':sha256_file(ppath),'all_union_regions_sha256':sha256_file(arpath)}
    write_json(root/'outputs/STAGE04A_FINAL_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
