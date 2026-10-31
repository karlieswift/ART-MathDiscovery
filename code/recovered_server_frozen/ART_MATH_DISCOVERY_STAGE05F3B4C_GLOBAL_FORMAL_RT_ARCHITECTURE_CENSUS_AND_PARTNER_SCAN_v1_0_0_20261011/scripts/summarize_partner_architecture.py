#!/usr/bin/env python3
import csv, json, hashlib, statistics, sys
from pathlib import Path
from collections import defaultdict, Counter

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()
def boolv(x): return str(x).lower()=='true'

def main():
    if len(sys.argv)!=8:
        raise SystemExit('usage: summarize_partner_architecture.py LOCUS_CENSUS ARRAYS NEAREST_ORF CLUSTER_TSV FAMILY_REPEAT_SUMMARY OUTDIR LOCK_JSON')
    locp,arrp,nearp,clustp,famrep,outdir,lockp=sys.argv[1:]
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True); lock=json.load(open(lockp)); assert lock['gate_pass'] is True
    loc=list(csv.DictReader(open(locp,encoding='utf-8'),delimiter='\t'))
    arr=list(csv.DictReader(open(arrp,encoding='utf-8'),delimiter='\t'))
    near=list(csv.DictReader(open(nearp,encoding='utf-8'),delimiter='\t'))
    fsum=list(csv.DictReader(open(famrep,encoding='utf-8'),delimiter='\t'))
    if len(loc)!=10598 or len(fsum)!=4912: raise SystemExit('FAIL input dimensions')

    member_to_rep={}
    if Path(clustp).exists() and Path(clustp).stat().st_size:
        with open(clustp,encoding='utf-8') as f:
            for line in f:
                if not line.strip(): continue
                x=line.rstrip('\n').split('\t')
                if len(x)>=2: member_to_rep[x[1]]=x[0]
    for n in near:
        n['partner_cluster_id']=member_to_rep.get(n['orf_id'],'') if n['orf_id'] else ''

    loc_by={r['candidate_id']:r for r in loc}
    strict_by=defaultdict(list); broad_by=defaultdict(list); up_by=defaultdict(list)
    for a in arr:
        if a['region']=='UPSTREAM': up_by[a['candidate_id']].append(a)
        if boolv(a['strict_published_art_geometry']): strict_by[a['candidate_id']].append(a)
        if boolv(a['broad_art_like_geometry']): broad_by[a['candidate_id']].append(a)
    near_by={r['candidate_id']:r for r in near}

    locus_arch=[]
    for cid,r in loc_by.items():
        n=near_by.get(cid,{})
        locus_arch.append({
          'family_rep_id':r['family_rep_id'],'candidate_id':cid,'accession':r['accession'],'placement_status':r['placement_status'],
          'has_upstream_repeat_array':bool(up_by.get(cid)),'has_strict_published_art_geometry':bool(strict_by.get(cid)),'has_broad_art_like_geometry':bool(broad_by.get(cid)),
          'nearest_downstream_orf_id':n.get('orf_id',''),'nearest_downstream_partner_cluster_id':n.get('partner_cluster_id',''),
          'nearest_orf_orientation_relative_to_rt':n.get('orientation_relative_to_rt',''),'nearest_orf_gap_after_rt_stop_bp':n.get('gap_after_rt_stop_bp',''),
          'nearest_orf_aa_length':n.get('aa_length',''),'rt_terminal_stop_present':n.get('rt_terminal_stop_present',''),
          'downstream_partner_scan_available_bp':n.get('partner_scan_available_bp',''),
        })
    lap=out/'F3B4C2_RT_REPEAT_PARTNER_LOCUS_ARCHITECTURE_10598.tsv'; lcols=list(locus_arch[0])
    with lap.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=lcols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(locus_arch)

    byfam=defaultdict(list)
    for r in locus_arch: byfam[r['family_rep_id']].append(r)
    famorder=[r['family_rep_id'] for r in fsum]
    rows=[]
    for fam in famorder:
        lr=byfam[fam]
        withp=[x for x in lr if x['nearest_downstream_orf_id']]
        def cluster_stats(sub):
            d=defaultdict(list)
            for x in sub:
                if x['nearest_downstream_partner_cluster_id']: d[x['nearest_downstream_partner_cluster_id']].append(x)
            if not d: return ('',0,0,0.0,'')
            ranked=sorted(d.items(),key=lambda kv:(-len({x['accession'] for x in kv[1]}),-len(kv[1]),kv[0]))
            cid,xs=ranked[0]; acc=len({x['accession'] for x in xs}); n=len(xs); denom=len(sub); frac=n/denom if denom else 0.0
            gaps=[int(x['nearest_orf_gap_after_rt_stop_bp']) for x in xs if str(x['nearest_orf_gap_after_rt_stop_bp'])!='']
            med=statistics.median(gaps) if gaps else ''
            return cid,n,acc,frac,med
        allc=cluster_stats(withp)
        up=[x for x in lr if boolv(x['has_upstream_repeat_array']) and x['nearest_downstream_orf_id']]
        strict=[x for x in lr if boolv(x['has_strict_published_art_geometry']) and x['nearest_downstream_orf_id']]
        broad=[x for x in lr if boolv(x['has_broad_art_like_geometry']) and x['nearest_downstream_orf_id']]
        uc=cluster_stats(up); sc=cluster_stats(strict); bc=cluster_stats(broad)
        rows.append({
          'family_rep_id':fam,'formal_locus_n':len(lr),'accession_n':len({x['accession'] for x in lr}),
          'locus_n_with_nearest_downstream_orf':len(withp),'accession_n_with_nearest_downstream_orf':len({x['accession'] for x in withp}),
          'same_strand_nearest_orf_locus_n':sum(x['nearest_orf_orientation_relative_to_rt']=='SAME_STRAND' for x in withp),
          'dominant_partner_cluster_id':allc[0],'dominant_partner_cluster_instance_n':allc[1],'dominant_partner_cluster_accession_n':allc[2],'dominant_partner_cluster_fraction_of_partner_loci':allc[3],'dominant_partner_cluster_median_gap_bp':allc[4],
          'upstream_array_locus_n':sum(boolv(x['has_upstream_repeat_array']) for x in lr),'upstream_array_locus_n_with_partner':len(up),
          'upstream_array_dominant_partner_cluster_id':uc[0],'upstream_array_dominant_partner_cluster_instance_n':uc[1],'upstream_array_dominant_partner_cluster_accession_n':uc[2],'upstream_array_dominant_partner_cluster_median_gap_bp':uc[4],
          'strict_array_locus_n':sum(boolv(x['has_strict_published_art_geometry']) for x in lr),'strict_array_locus_n_with_partner':len(strict),
          'strict_array_dominant_partner_cluster_id':sc[0],'strict_array_dominant_partner_cluster_instance_n':sc[1],'strict_array_dominant_partner_cluster_accession_n':sc[2],'strict_array_dominant_partner_cluster_median_gap_bp':sc[4],
          'strict_array_plus_recurrent_partner_cluster':bool(sc[2]>=2),
          'broad_array_locus_n':sum(boolv(x['has_broad_art_like_geometry']) for x in lr),'broad_array_locus_n_with_partner':len(broad),
          'broad_array_dominant_partner_cluster_id':bc[0],'broad_array_dominant_partner_cluster_instance_n':bc[1],'broad_array_dominant_partner_cluster_accession_n':bc[2],'broad_array_dominant_partner_cluster_median_gap_bp':bc[4],
          'broad_array_plus_recurrent_partner_cluster':bool(bc[2]>=2),
          'partner_function_confirmed':False,'dedicated_partner_gene_confirmed':False,
        })
    fp=out/'F3B4C2_RT_FAMILY_REPEAT_PARTNER_ARCHITECTURE_SUMMARY_4912.tsv'; fcols=list(rows[0])
    with fp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fcols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)

    cand=[r for r in rows if r['strict_array_locus_n'] or r['broad_array_locus_n'] or r['upstream_array_locus_n']]
    cand.sort(key=lambda r:(-int(r['strict_array_dominant_partner_cluster_accession_n']),-int(r['strict_array_locus_n']),-int(r['broad_array_dominant_partner_cluster_accession_n']),-int(r['upstream_array_dominant_partner_cluster_accession_n']),r['family_rep_id']))
    cp=out/'F3B4C2_ARRAY_ASSOCIATED_RT_PARTNER_CANDIDATE_FAMILIES.tsv'
    with cp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fcols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(cand)

    gate={
      'stage':'Stage05F3B4C2','gate':'PASS_F3B4C2_GLOBAL_FORMAL_RT_DOWNSTREAM_PARTNER_ARCHITECTURE_SCAN_COMPLETE','gate_pass':True,
      'formal_rt_family_n':len(rows),'formal_rt_locus_n':len(locus_arch),
      'family_n_with_upstream_repeat_array':sum(r['upstream_array_locus_n']>0 for r in rows),
      'family_n_with_strict_repeat_geometry':sum(r['strict_array_locus_n']>0 for r in rows),
      'family_n_with_broad_repeat_geometry':sum(r['broad_array_locus_n']>0 for r in rows),
      'family_n_with_strict_repeat_and_recurrent_nearest_partner_cluster':sum(boolv(r['strict_array_plus_recurrent_partner_cluster']) for r in rows),
      'family_n_with_broad_repeat_and_recurrent_nearest_partner_cluster':sum(boolv(r['broad_array_plus_recurrent_partner_cluster']) for r in rows),
      'recurrent_partner_definition':'same operational partner-protein cluster observed as nearest downstream ORF in >=2 independent accessions among loci carrying the indicated repeat geometry',
      'partner_cluster_rule':'MMseqs2 operational clustering min-seq-id 0.30, coverage 0.70, cov-mode 0, cluster-mode 0, sensitivity 7.5',
      'nearest_orf_is_confirmed_partner_gene':False,'partner_function_decision_performed':False,'novelty_decision_performed':False,
      'locus_architecture_sha256':sha256_file(lap),'family_architecture_sha256':sha256_file(fp),'candidate_family_sha256':sha256_file(cp),
      'claim_boundary':'Computational partner-candidate architecture only. A recurrent nearest downstream ORF cluster strengthens ART-like architecture but does not by itself establish a dedicated partner gene, biochemical interaction, function, or novelty.'}
    gp=out/'STAGE05F3B4C2_GATE.json'; gp.write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2))
    print('\n=== ARRAY-ASSOCIATED FAMILY PARTNER SUMMARY ===')
    print('\t'.join(fcols))
    for r in cand[:100]: print('\t'.join(str(r[c]) for c in fcols))
    print('PASS_STAGE05F3B4C2_PARTNER_ARCHITECTURE_SCAN_COMPLETE')

if __name__=='__main__': main()
