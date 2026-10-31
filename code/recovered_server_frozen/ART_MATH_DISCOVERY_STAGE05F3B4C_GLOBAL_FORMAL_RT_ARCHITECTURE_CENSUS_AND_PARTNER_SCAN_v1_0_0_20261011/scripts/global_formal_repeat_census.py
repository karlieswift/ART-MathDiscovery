#!/usr/bin/env python3
import csv, json, hashlib, importlib.util, sys
from pathlib import Path
from collections import defaultdict


def loadmod(name, path):
    spec=importlib.util.spec_from_file_location(name, path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def boolv(x): return str(x).lower()=='true'

def main():
    if len(sys.argv)!=12:
        raise SystemExit('usage: global_formal_repeat_census.py R4_FAMILY_MAP R4_INSTANCE_MAP ORFDB ORFDB_INDEX GENOMEDB GENOMEDB_INDEX OUTDIR LOCK_JSON VENDOR_PLACEMENT VENDOR_REPEAT R4_FAMILY_SUMMARY')
    fmap_p,imap_p,orfdb,orfidx,gdb,gidx,outdir,lockp,vplace,vrepeat,fsump=sys.argv[1:]
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    lock=json.load(open(lockp)); assert lock['gate_pass'] is True
    placement=loadmod('vendor_place',vplace)
    rpt=loadmod('vendor_repeat',vrepeat)

    members=[]
    with open(fmap_p,encoding='utf-8') as f:
        rr=csv.DictReader(f,delimiter='\t')
        req={'candidate_id','family_rep_id','accession'}
        if not req.issubset(rr.fieldnames or []): raise SystemExit(f'FAIL family map missing {sorted(req-set(rr.fieldnames or []))}')
        for r in rr: members.append(r)
    if len(members)!=10598:
        raise SystemExit(f'FAIL expected 10598 formal loci got {len(members)}')
    famset={r['family_rep_id'] for r in members}
    if len(famset)!=4912: raise SystemExit(f'FAIL expected 4912 families got {len(famset)}')

    im={}
    wanted={r['candidate_id'] for r in members}
    with open(imap_p,encoding='utf-8') as f:
        rr=csv.DictReader(f,delimiter='\t')
        req={'instance_id','orfdb_key','parent_contig_key','rescued_accession','parent_contig_header'}
        if not req.issubset(rr.fieldnames or []): raise SystemExit(f'FAIL instance map missing {sorted(req-set(rr.fieldnames or []))}')
        for r in rr:
            if r['instance_id'] in wanted: im[r['instance_id']]=r
    miss=wanted-set(im)
    if miss: raise SystemExit(f'FAIL provenance missing n={len(miss)} examples={sorted(miss)[:10]}')

    orfkeys={im[r['candidate_id']]['orfdb_key'] for r in members}
    contigkeys={im[r['candidate_id']]['parent_contig_key'] for r in members}
    oi=placement.load_index(orfidx,orfkeys); gi=placement.load_index(gidx,contigkeys)

    arrays=[]; locus=[]
    with open(orfdb,'rb',buffering=0) as ofh, open(gdb,'rb',buffering=0) as gfh:
        for ix,m in enumerate(members,1):
            q=m['candidate_id']; p=im[q]
            pep=placement.fetch_db(ofh,*oi[p['orfdb_key']]).rstrip('*')
            seq=placement.fetch_db(gfh,*gi[p['parent_contig_key']])
            hits=placement.exact_hits(seq,pep)
            status='UNIQUE_EXACT' if len(hits)==1 else ('NO_EXACT' if not hits else 'MULTI_EXACT')
            base={
                'family_rep_id':m['family_rep_id'],'candidate_id':q,
                'accession':p['rescued_accession'],'parent_contig_key':p['parent_contig_key'],
                'parent_contig_header':p['parent_contig_header'],'parent_contig_length_bp':len(seq),
                'aa_length':len(pep),'exact_peptide_match_n':len(hits),'placement_status':status,
            }
            if len(hits)!=1:
                base.update({'strand':'','frame_zero_based':'','candidate_start0':'','candidate_end0':'',
                             'upstream_available_bp':0,'downstream_available_bp':0,
                             'any_generic_repeat_array':False,'upstream_repeat_array_n':0,'downstream_repeat_array_n':0,
                             'strict_published_art_array_n':0,'broad_art_like_array_n':0,
                             'upstream_complete_500bp':False,'upstream_complete_2000bp':False,'upstream_complete_5000bp':False,
                             'no_array_in_complete_5kb_upstream':False})
                locus.append(base); continue
            strand,frame,s,e=hits[0]
            if strand=='+':
                up=seq[max(0,s-rpt.SCAN_FLANK):s]
                dn=seq[e:min(len(seq),e+rpt.SCAN_FLANK)]
            else:
                up=placement.revcomp(seq[e:min(len(seq),e+rpt.SCAN_FLANK)])
                dn=placement.revcomp(seq[max(0,s-rpt.SCAN_FLANK):s])
            up_arr=rpt.scan_region(up,'UPSTREAM',-len(up))
            dn_arr=rpt.scan_region(dn,'DOWNSTREAM',len(pep)*3)
            arr=up_arr+dn_arr
            for j,a in enumerate(arr,1):
                z={'family_rep_id':m['family_rep_id'],'candidate_id':q,'accession':p['rescued_accession'],
                   'parent_contig_key':p['parent_contig_key'],'strand':strand,
                   'upstream_available_bp':len(up),'downstream_available_bp':len(dn),'array_rank_within_locus':j}
                z.update(a); arrays.append(z)
            base.update({
                'strand':strand,'frame_zero_based':frame,'candidate_start0':s,'candidate_end0':e,
                'upstream_available_bp':len(up),'downstream_available_bp':len(dn),
                'any_generic_repeat_array':bool(arr),'upstream_repeat_array_n':len(up_arr),'downstream_repeat_array_n':len(dn_arr),
                'strict_published_art_array_n':sum(boolv(a['strict_published_art_geometry']) for a in arr),
                'broad_art_like_array_n':sum(boolv(a['broad_art_like_geometry']) for a in arr),
                'upstream_complete_500bp':len(up)>=500,'upstream_complete_2000bp':len(up)>=2000,'upstream_complete_5000bp':len(up)>=5000,
                'no_array_in_complete_5kb_upstream':bool(len(up)>=5000 and not up_arr),
            })
            locus.append(base)
            if ix%250==0:
                print(f'[SCAN] {ix}/{len(members)} arrays={len(arrays)}',flush=True)

    locp=out/'F3B4C1_GLOBAL_FORMAL_RT_LOCUS_REPEAT_CENSUS_10598.tsv'
    loccols=list(locus[0])
    with locp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=loccols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(locus)

    arrp=out/'F3B4C1_ALL_DETECTED_REPEAT_ARRAYS.tsv'
    if arrays:
        arrcols=list(arrays[0])
    else:
        arrcols=['family_rep_id','candidate_id','accession','parent_contig_key','strand','upstream_available_bp','downstream_available_bp','array_rank_within_locus','region','seed_k','motif','motif_len','copy_n','local_repeat_starts','candidate_relative_repeat_starts','start_to_start_gaps','spacer_lengths','median_spacer_bp','min_spacer_bp','max_spacer_bp','spacer_cv','array_span_bp','strict_published_art_geometry','broad_art_like_geometry']
    with arrp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=arrcols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(arrays)

    famorder=[]
    with open(fsump,encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'): famorder.append(r['family_rep_id'])
    if len(famorder)!=4912 or set(famorder)!=famset:
        famorder=sorted(famset)

    byfam=defaultdict(list)
    for r in locus: byfam[r['family_rep_id']].append(r)
    arrby=defaultdict(list)
    for a in arrays: arrby[a['family_rep_id']].append(a)
    fs=[]
    for fam in famorder:
        lr=byfam[fam]; ar=arrby[fam]
        strict=[a for a in ar if boolv(a['strict_published_art_geometry'])]
        broad=[a for a in ar if boolv(a['broad_art_like_geometry'])]
        up=[a for a in ar if a['region']=='UPSTREAM']
        fs.append({
            'family_rep_id':fam,'formal_locus_n':len(lr),'accession_n':len({x['accession'] for x in lr}),
            'unique_exact_locus_n':sum(x['placement_status']=='UNIQUE_EXACT' for x in lr),
            'instance_n_with_any_repeat_array':sum(boolv(x['any_generic_repeat_array']) for x in lr),
            'accession_n_with_any_repeat_array':len({x['accession'] for x in lr if boolv(x['any_generic_repeat_array'])}),
            'instance_n_with_upstream_repeat_array':len({a['candidate_id'] for a in up}),
            'accession_n_with_upstream_repeat_array':len({a['accession'] for a in up}),
            'strict_published_art_array_n':len(strict),'strict_published_art_instance_n':len({a['candidate_id'] for a in strict}),
            'strict_published_art_accession_n':len({a['accession'] for a in strict}),
            'broad_art_like_array_n':len(broad),'broad_art_like_instance_n':len({a['candidate_id'] for a in broad}),
            'broad_art_like_accession_n':len({a['accession'] for a in broad}),
            'complete_5kb_upstream_instance_n':sum(boolv(x['upstream_complete_5000bp']) for x in lr),
            'complete_5kb_upstream_no_array_instance_n':sum(boolv(x['no_array_in_complete_5kb_upstream']) for x in lr),
            'architecture_absence_claim_authorized':False,
        })
    fsp=out/'F3B4C1_GLOBAL_FORMAL_RT_FAMILY_REPEAT_SUMMARY_4912.tsv'; fscols=list(fs[0])
    with fsp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fscols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(fs)

    strict_rows=[]
    strict_candidates={a['candidate_id'] for a in arrays if boolv(a['strict_published_art_geometry'])}
    broad_candidates={a['candidate_id'] for a in arrays if boolv(a['broad_art_like_geometry'])}
    up_candidates={a['candidate_id'] for a in arrays if a['region']=='UPSTREAM'}
    for r in locus:
        if r['candidate_id'] in strict_candidates or r['candidate_id'] in broad_candidates or r['candidate_id'] in up_candidates:
            z=dict(r); z['has_upstream_array']=r['candidate_id'] in up_candidates; z['has_strict_published_art_geometry']=r['candidate_id'] in strict_candidates; z['has_broad_art_like_geometry']=r['candidate_id'] in broad_candidates; strict_rows.append(z)
    candp=out/'F3B4C1_UPSTREAM_REPEAT_AND_ART_GEOMETRY_LOCUS_CANDIDATES.tsv'
    ccols=list(strict_rows[0]) if strict_rows else loccols+['has_upstream_array','has_strict_published_art_geometry','has_broad_art_like_geometry']
    with candp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=ccols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(strict_rows)

    gate={
      'stage':'Stage05F3B4C1','gate':'PASS_F3B4C1_GLOBAL_10598_FORMAL_RT_REPEAT_ARCHITECTURE_CENSUS_COMPLETE','gate_pass':True,
      'formal_rt_locus_n':len(locus),'formal_rt_family_n':len(fs),
      'unique_exact_placement_n':sum(x['placement_status']=='UNIQUE_EXACT' for x in locus),
      'no_exact_placement_n':sum(x['placement_status']=='NO_EXACT' for x in locus),
      'multi_exact_placement_n':sum(x['placement_status']=='MULTI_EXACT' for x in locus),
      'detected_repeat_array_n':len(arrays),
      'family_n_with_any_repeat_array':sum(x['instance_n_with_any_repeat_array']>0 for x in fs),
      'family_n_with_upstream_repeat_array':sum(x['instance_n_with_upstream_repeat_array']>0 for x in fs),
      'family_n_with_strict_published_art_geometry':sum(x['strict_published_art_array_n']>0 for x in fs),
      'family_n_with_broad_art_like_geometry':sum(x['broad_art_like_array_n']>0 for x in fs),
      'strict_published_art_instance_n':len(strict_candidates),'broad_art_like_instance_n':len(broad_candidates),
      'repeat_detector_parameters_changed':False,'repeat_detector_role':'byte-identical F3B4B1 detector applied exhaustively to all 10,598 formal RT loci',
      'architecture_absence_claim_authorized':False,'partner_analysis_performed':False,'novelty_decision_performed':False,
      'locus_census_sha256':sha256_file(locp),'all_arrays_sha256':sha256_file(arrp),'family_summary_sha256':sha256_file(fsp),'partner_handoff_sha256':sha256_file(candp),
      'claim_boundary':'Exhaustive repeat-architecture census across the frozen 10,598 formal RVT-supported loci. ART geometry is descriptive annotation only; partner gene and novelty are separate downstream evidence.'}
    gp=out/'STAGE05F3B4C1_GATE.json'; gp.write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2))
    print('\n=== STRICT/BROAD FAMILY TOPLINES ===')
    interesting=[x for x in fs if x['strict_published_art_array_n'] or x['broad_art_like_array_n']]
    interesting.sort(key=lambda x:(-x['strict_published_art_accession_n'],-x['strict_published_art_instance_n'],-x['broad_art_like_accession_n'],-x['instance_n_with_upstream_repeat_array'],x['family_rep_id']))
    for x in interesting[:100]: print('\t'.join(str(x[c]) for c in fscols))
    print('PASS_STAGE05F3B4C1_GLOBAL_REPEAT_CENSUS_COMPLETE')

if __name__=='__main__': main()
