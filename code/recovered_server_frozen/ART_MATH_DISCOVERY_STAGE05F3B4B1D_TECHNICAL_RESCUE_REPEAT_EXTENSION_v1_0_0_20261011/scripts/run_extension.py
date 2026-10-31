#!/usr/bin/env python3
import csv, json, hashlib, importlib.util, sys
from pathlib import Path
from collections import defaultdict


def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def loadmod(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def boolv(x): return str(x).lower()=='true'

def write_tsv(path, rows, cols=None):
    path=Path(path)
    if cols is None: cols=list(rows[0]) if rows else []
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def main():
    if len(sys.argv)!=15:
        raise SystemExit('usage: run_extension.py ER1_SUMMARY FM IM ORFDB ORFIDX GDB GIDX B4B1_FAMSUM B4B1C_WEAK_FAMSUM OUTDIR LOCK DETECTOR_PY PLACEMENT_PY EXPECTED_TARGET_N')
    er1p,fmp,imp,orfdb,orfidx,gdb,gidx,b4b1sum,b4b1cweak,outdir,lockp,detpy,placepy,nexp=sys.argv[1:]
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    lock=json.load(open(lockp)); assert lock['gate_pass'] is True
    detector=loadmod('frozen_b4b1_repeat',detpy)
    placement=loadmod('frozen_b4b0_place',placepy)

    er1=list(csv.DictReader(open(er1p,encoding='utf-8'),delimiter='\t'))
    if len(er1)!=32: raise SystemExit(f'FAIL ER1 rows={len(er1)}')
    targets=[r for r in er1 if r.get('b3er1_resolution_source')=='FOLDSEEK_ONLY_TECHNICAL_REPAIR']
    targets=sorted(targets,key=lambda r:int(r['shortlist_rank']))
    if len(targets)!=int(nexp): raise SystemExit(f'FAIL repaired target n={len(targets)} expected={nexp}')
    finaltiers=[r['b3er1_final_best_tier'] for r in targets]
    if sum(t in {'STRONG','SUGGESTIVE'} for t in finaltiers)!=2 or sum(t=='WEAK' for t in finaltiers)!=3:
        raise SystemExit(f'FAIL expected repaired final tiers 2 positive + 3 weak, got {finaltiers}')
    fam_order=[r['family_rep_id'] for r in targets]
    famset=set(fam_order); tier={r['family_rep_id']:r['b3er1_final_best_tier'] for r in targets}; pri={r['family_rep_id']:r['candidate_id'] for r in targets}; rank={r['family_rep_id']:r['shortlist_rank'] for r in targets}

    members=[]
    with open(fmp,encoding='utf-8') as f:
        rr=csv.DictReader(f,delimiter='\t')
        for r in rr:
            if r['family_rep_id'] in famset: members.append(r)
    byfam=defaultdict(list)
    for r in members: byfam[r['family_rep_id']].append(r)
    miss=famset-set(byfam)
    if miss: raise SystemExit(f'FAIL no members for repaired families {sorted(miss)}')

    im={}; wanted={r['candidate_id'] for r in members}
    with open(imp,encoding='utf-8') as f:
        rr=csv.DictReader(f,delimiter='\t')
        for r in rr:
            if r['instance_id'] in wanted: im[r['instance_id']]=r
    missing=wanted-set(im)
    if missing: raise SystemExit(f'FAIL missing provenance n={len(missing)} examples={sorted(missing)[:10]}')

    orfkeys={im[r['candidate_id']]['orfdb_key'] for r in members}; contigkeys={im[r['candidate_id']]['parent_contig_key'] for r in members}
    oi=placement.load_index(orfidx,orfkeys); gi=placement.load_index(gidx,contigkeys)
    peps={}; seqs={}
    with open(orfdb,'rb',buffering=0) as fh:
        for k,(o,l) in oi.items(): peps[k]=placement.fetch_db(fh,o,l).rstrip('*')
    with open(gdb,'rb',buffering=0) as fh:
        for k,(o,l) in gi.items(): seqs[k]=placement.fetch_db(fh,o,l)

    loci=[]; arrays=[]
    for ix,m in enumerate(members,1):
        q=m['candidate_id']; p=im[q]; pep=peps[p['orfdb_key']]; seq=seqs[p['parent_contig_key']]
        hits=placement.exact_hits(seq,pep)
        if len(hits)!=1: raise SystemExit(f'FAIL repaired-family exact placement {q} n={len(hits)}')
        strand,frame,s,e=hits[0]
        if strand=='+':
            up=seq[max(0,s-detector.SCAN_FLANK):s]; dn=seq[e:min(len(seq),e+detector.SCAN_FLANK)]
        else:
            up=detector.revcomp(seq[e:min(len(seq),e+detector.SCAN_FLANK)]); dn=detector.revcomp(seq[max(0,s-detector.SCAN_FLANK):s])
        up_av=len(up); dn_av=len(dn)
        up_arr=detector.scan_region(up,'UPSTREAM',-up_av); dn_arr=detector.scan_region(dn,'DOWNSTREAM',len(pep)*3); arr=up_arr+dn_arr
        for j,a in enumerate(arr,1):
            z={'family_rep_id':m['family_rep_id'],'structural_tier':tier[m['family_rep_id']],'shortlist_rank':rank[m['family_rep_id']],
               'priority_candidate_id':pri[m['family_rep_id']],'candidate_id':q,'accession':p['rescued_accession'],'parent_contig_key':p['parent_contig_key'],
               'strand':strand,'upstream_available_bp':up_av,'downstream_available_bp':dn_av,'array_rank_within_locus':j}
            z.update(a); arrays.append(z)
        loci.append({'family_rep_id':m['family_rep_id'],'structural_tier':tier[m['family_rep_id']],'shortlist_rank':rank[m['family_rep_id']],
                     'priority_candidate_id':pri[m['family_rep_id']],'candidate_id':q,'accession':p['rescued_accession'],'parent_contig_key':p['parent_contig_key'],
                     'parent_contig_length_bp':len(seq),'aa_length':len(pep),'strand':strand,'frame_zero_based':frame,'candidate_start0':s,'candidate_end0':e,
                     'left_available_bp':s,'right_available_bp':len(seq)-e,'min_flank_bp':min(s,len(seq)-e),'total_flank_bp':s+len(seq)-e,
                     'upstream_available_bp':up_av,'downstream_available_bp':dn_av,'any_generic_repeat_array':bool(arr),'upstream_repeat_array_n':len(up_arr),
                     'downstream_repeat_array_n':len(dn_arr),'strict_published_art_array_n':sum(bool(a['strict_published_art_geometry']) for a in arr),
                     'broad_art_like_array_n':sum(bool(a['broad_art_like_geometry']) for a in arr),
                     'no_array_in_complete_5kb_upstream':bool(up_av>=detector.ABSENCE_UPSTREAM_COMPLETE and not up_arr)})
        if ix%100==0: print(f'[SCAN] {ix}/{len(members)} arrays={len(arrays)}',flush=True)

    locp=out/'F3B4B1D_REPAIRED5_ALL_MEMBER_LOCUS_REPEAT_SUMMARY.tsv'; arrp=out/'F3B4B1D_REPAIRED5_ALL_DETECTED_REPEAT_ARRAYS.tsv'
    write_tsv(locp,loci)
    arrcols=list(arrays[0]) if arrays else ['family_rep_id','structural_tier','shortlist_rank','priority_candidate_id','candidate_id','accession','parent_contig_key','strand','upstream_available_bp','downstream_available_bp','array_rank_within_locus','region','seed_k','motif','motif_len','copy_n','local_repeat_starts','candidate_relative_repeat_starts','start_to_start_gaps','spacer_lengths','median_spacer_bp','min_spacer_bp','max_spacer_bp','spacer_cv','array_span_bp','strict_published_art_geometry','broad_art_like_geometry']
    write_tsv(arrp,arrays,arrcols)

    famsum=[]
    for fam in fam_order:
        lr=[x for x in loci if x['family_rep_id']==fam]; ar=[x for x in arrays if x['family_rep_id']==fam]
        strict=[x for x in ar if boolv(x['strict_published_art_geometry'])]; broad=[x for x in ar if boolv(x['broad_art_like_geometry'])]; up=[x for x in ar if x['region']=='UPSTREAM']
        famsum.append({'family_rep_id':fam,'structural_tier':tier[fam],'shortlist_rank':rank[fam],'priority_candidate_id':pri[fam],
          'family_instance_n':len(lr),'family_accession_n':len({x['accession'] for x in lr}),
          'instance_n_with_any_repeat_array':sum(bool(x['any_generic_repeat_array']) for x in lr),'accession_n_with_any_repeat_array':len({x['accession'] for x in lr if x['any_generic_repeat_array']}),
          'instance_n_with_upstream_repeat_array':len({x['candidate_id'] for x in up}),'accession_n_with_upstream_repeat_array':len({x['accession'] for x in up}),
          'strict_published_art_array_n':len(strict),'strict_published_art_instance_n':len({x['candidate_id'] for x in strict}),'strict_published_art_accession_n':len({x['accession'] for x in strict}),
          'broad_art_like_array_n':len(broad),'broad_art_like_instance_n':len({x['candidate_id'] for x in broad}),'broad_art_like_accession_n':len({x['accession'] for x in broad}),
          'complete_5kb_upstream_instance_n':sum(int(x['upstream_available_bp'])>=5000 for x in lr),
          'complete_5kb_upstream_no_array_instance_n':sum(bool(x['no_array_in_complete_5kb_upstream']) for x in lr),
          'max_min_flank_bp':max(int(x['min_flank_bp']) for x in lr) if lr else 0})
    fsp=out/'F3B4B1D_REPAIRED5_FAMILY_REPEAT_SUMMARY.tsv'; write_tsv(fsp,famsum)

    # Integrate original 16 positives + original 9 WEAK + repaired extension 2 positive + 3 WEAK.
    sv=list(csv.DictReader(open(b4b1sum,encoding='utf-8'),delimiter='\t')); wk=list(csv.DictReader(open(b4b1cweak,encoding='utf-8'),delimiter='\t'))
    if len(sv)!=16 or len(wk)!=9: raise SystemExit(f'FAIL base summaries survivor={len(sv)} weak={len(wk)}')
    def canon(r,source,default_tier=None):
        t=r.get('structural_tier') or r.get('survivor_tier') or default_tier or ''
        z=dict(r); z['structural_tier']=t; z['architecture_source']=source; return z
    integrated=[canon(r,'B4B1_ORIGINAL16') for r in sv]+[canon(r,'B4B1C_ORIGINAL_WEAK9','WEAK') for r in wk]+[canon(r,'B4B1D_TECHNICAL_RESCUE5') for r in famsum]
    if len(integrated)!=30: raise SystemExit('FAIL integrated resolved family n')
    if len({r['family_rep_id'] for r in integrated})!=30: raise SystemExit('FAIL duplicate integrated family IDs')
    pos=[r for r in integrated if r['structural_tier'] in {'STRONG','SUGGESTIVE'}]; weak=[r for r in integrated if r['structural_tier']=='WEAK']
    if len(pos)!=18 or len(weak)!=12: raise SystemExit(f'FAIL final groups positive={len(pos)} weak={len(weak)}')
    intp=out/'F3B4B1D_FINAL_RESOLVED30_FAMILY_REPEAT_SUMMARY.tsv'
    keys=['family_rep_id','structural_tier','architecture_source','family_instance_n','family_accession_n','instance_n_with_any_repeat_array','accession_n_with_any_repeat_array','instance_n_with_upstream_repeat_array','accession_n_with_upstream_repeat_array','strict_published_art_array_n','strict_published_art_instance_n','strict_published_art_accession_n','broad_art_like_array_n','broad_art_like_instance_n','broad_art_like_accession_n']
    norm=[]
    for r in integrated: norm.append({k:r.get(k,'') for k in keys})
    write_tsv(intp,norm,keys)

    metrics=[]
    for metric in ['instance_n_with_any_repeat_array','instance_n_with_upstream_repeat_array','strict_published_art_instance_n','broad_art_like_instance_n']:
        pn=sum(int(r[metric]) for r in pos); pd=sum(int(r['family_instance_n']) for r in pos); wn=sum(int(r[metric]) for r in weak); wd=sum(int(r['family_instance_n']) for r in weak)
        metrics.append({'metric':metric,'strong_suggestive_positive_n':pn,'strong_suggestive_denominator_n':pd,'strong_suggestive_rate':pn/pd if pd else 0,
                        'weak_positive_n':wn,'weak_denominator_n':wd,'weak_rate':wn/wd if wd else 0,'rate_ratio':(pn/pd)/(wn/wd) if pd and wd and wn else ('INF' if pn else 0)})
    for metric in ['instance_n_with_any_repeat_array','instance_n_with_upstream_repeat_array','strict_published_art_array_n','broad_art_like_array_n']:
        pn=sum(int(r[metric])>0 for r in pos); wn=sum(int(r[metric])>0 for r in weak)
        metrics.append({'metric':'family_any__'+metric,'strong_suggestive_positive_n':pn,'strong_suggestive_denominator_n':18,'strong_suggestive_rate':pn/18,
                        'weak_positive_n':wn,'weak_denominator_n':12,'weak_rate':wn/12,'rate_ratio':(pn/18)/(wn/12) if wn else ('INF' if pn else 0)})
    cp=out/'F3B4B1D_FINAL_RESOLVED30_DESCRIPTIVE_COMPARISON.tsv'; write_tsv(cp,metrics)

    gate={'stage':'Stage05F3B4B1D','gate':'PASS_F3B4B1D_TECHNICAL_RESCUE_REPEAT_EXTENSION_COMPLETE','gate_pass':True,
          'analysis_role':'post-hoc technical-rescue extension using unchanged B4B1 repeat detector; original B4B1/B4B1C outputs preserved',
          'repaired_family_n':5,'repaired_strong_or_suggestive_family_n':2,'repaired_weak_family_n':3,'repaired_member_locus_n':len(loci),'repaired_detected_array_n':len(arrays),
          'repaired_family_n_with_any_repeat_array':sum(int(r['instance_n_with_any_repeat_array'])>0 for r in famsum),
          'repaired_family_n_with_upstream_repeat_array':sum(int(r['instance_n_with_upstream_repeat_array'])>0 for r in famsum),
          'repaired_family_n_with_strict_geometry':sum(int(r['strict_published_art_array_n'])>0 for r in famsum),
          'repaired_family_n_with_broad_geometry':sum(int(r['broad_art_like_array_n'])>0 for r in famsum),
          'integrated_resolved_family_n':30,'integrated_strong_or_suggestive_family_n':18,'integrated_weak_family_n':12,'resource_limited_unresolved_family_n':2,
          'repeat_detector_parameters_changed':False,'original_b4b1_membership_changed':False,'original_b4b1c_control_changed':False,'unresolved_treated_as_negative':False,'novelty_decision_performed':False,
          'repaired_family_summary_sha256':sha256_file(fsp),'integrated_family_summary_sha256':sha256_file(intp),'integrated_comparison_sha256':sha256_file(cp),
          'claim_boundary':'Technical-rescue architecture extension only. Final 30-family repeat comparison is descriptive/post-discovery; two resource-limited families remain censored.'}
    gp=out/'STAGE05F3B4B1D_GATE.json'; gp.write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2))
    print('\n=== REPAIRED FIVE FAMILY ARCHITECTURE SUMMARY ===')
    print('\t'.join(famsum[0].keys()))
    for r in famsum: print('\t'.join(str(r[k]) for k in famsum[0].keys()))
    print('\n=== FINAL RESOLVED30 DESCRIPTIVE COMPARISON ===')
    print('\t'.join(metrics[0].keys()))
    for r in metrics: print('\t'.join(str(r[k]) for k in metrics[0].keys()))
    print('PASS_STAGE05F3B4B1D_TECHNICAL_RESCUE_REPEAT_EXTENSION_COMPLETE')

if __name__=='__main__': main()
