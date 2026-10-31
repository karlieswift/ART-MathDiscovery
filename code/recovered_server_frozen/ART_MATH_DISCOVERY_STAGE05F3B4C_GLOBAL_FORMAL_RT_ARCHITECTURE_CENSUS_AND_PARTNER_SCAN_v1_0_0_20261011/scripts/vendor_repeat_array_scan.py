#!/usr/bin/env python3
import csv, json, hashlib, math, re, statistics, sys
from pathlib import Path
from collections import defaultdict, Counter

RCMAP=str.maketrans('ACGTRYMKSWBDHVNacgtrymkswbdhvn','TGCAYRKMSWVHDBNtgcayrkmswvhdbn')
SEED_K=(8,10,12,14,16)
MIN_COPIES=3
GAP_MIN=60
GAP_MAX=360
MAX_MOTIF=80
SCAN_FLANK=10000
ABSENCE_UPSTREAM_COMPLETE=5000

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def clean_entry(b): return b.rstrip(b'\x00\r\n')
def revcomp(s): return s.translate(RCMAP)[::-1]
def fetch_db(dbfh,off,ln):
    dbfh.seek(off); b=dbfh.read(ln)
    if len(b)!=ln: raise SystemExit(f'FAIL short DB read off={off} got={len(b)} expected={ln}')
    return re.sub(r'\s+','',clean_entry(b).decode('ascii','strict')).upper()
def load_index(path,wanted):
    out={}
    with open(path,encoding='ascii') as f:
        for line in f:
            if not line.strip(): continue
            x=line.split()
            if len(x)<3: raise SystemExit(f'FAIL malformed index {line[:120]!r}')
            if x[0] in wanted: out[x[0]]=(int(x[1]),int(x[2]))
    miss=wanted-set(out)
    if miss: raise SystemExit(f'FAIL missing DB keys n={len(miss)} examples={sorted(miss)[:10]}')
    return out

def entropy(s):
    n=len(s)
    if not n: return 0.0
    c=Counter(s)
    return -sum((v/n)*math.log2(v/n) for v in c.values())
def low_complexity(s):
    if not s: return True
    m=max(Counter(s).values())/len(s)
    return entropy(s)<1.35 or m>0.72

def chains_from_positions(pos):
    pos=sorted(pos)
    chains=[]; cur=[]
    for p in pos:
        if not cur:
            cur=[p]
        elif GAP_MIN <= p-cur[-1] <= GAP_MAX:
            cur.append(p)
        else:
            if len(cur)>=MIN_COPIES: chains.append(cur)
            cur=[p]
    if len(cur)>=MIN_COPIES: chains.append(cur)
    return chains

def extend_exact_core(seq, positions, k):
    # Extend the exact common seed left/right across all copies; max motif length frozen.
    left=0; right=k
    while left>-MAX_MOTIF and -left + right < MAX_MOTIF:
        idx=[p+left-1 for p in positions]
        if min(idx)<0: break
        b={seq[i] for i in idx}
        if len(b)!=1: break
        left-=1
    while -left + right < MAX_MOTIF:
        idx=[p+right for p in positions]
        if max(idx)>=len(seq): break
        b={seq[i] for i in idx}
        if len(b)!=1: break
        right+=1
    starts=[p+left for p in positions]
    motif=seq[starts[0]:starts[0]+(-left+right)]
    return motif,starts

def scan_region(seq, region_label, offset_rel_candidate):
    # Returns arrays in candidate-oriented sequence; rel coordinates use RT start=0, upstream negative.
    raw=[]
    for k in SEED_K:
        d=defaultdict(list)
        for i in range(0,len(seq)-k+1):
            mer=seq[i:i+k]
            if 'N' in mer or low_complexity(mer): continue
            d[mer].append(i)
        for mer,pos in d.items():
            if len(pos)<MIN_COPIES: continue
            for chain in chains_from_positions(pos):
                motif,starts=extend_exact_core(seq,chain,k)
                if len(motif)<8 or low_complexity(motif): continue
                gaps=[starts[i+1]-starts[i] for i in range(len(starts)-1)]
                spacers=[g-len(motif) for g in gaps]
                if any(x<0 for x in spacers): continue
                mean=statistics.mean(spacers) if spacers else 0.0
                med=statistics.median(spacers) if spacers else 0.0
                sd=statistics.pstdev(spacers) if len(spacers)>1 else 0.0
                cv=(sd/mean) if mean>0 else 0.0
                rel=[offset_rel_candidate+x for x in starts]
                strict_art=(region_label=='UPSTREAM' and 3<=len(starts)<=21 and all(120<=x<=220 for x in spacers))
                broad_art=(region_label=='UPSTREAM' and 3<=len(starts)<=30 and 100<=med<=240 and (sum(80<=x<=280 for x in spacers)/len(spacers) if spacers else 0)>=0.75)
                raw.append({
                    'region':region_label,'seed_k':k,'motif':motif,'motif_len':len(motif),'copy_n':len(starts),
                    'local_repeat_starts':','.join(map(str,starts)),'candidate_relative_repeat_starts':','.join(map(str,rel)),
                    'start_to_start_gaps':','.join(map(str,gaps)),'spacer_lengths':','.join(map(str,spacers)),
                    'median_spacer_bp':med,'min_spacer_bp':min(spacers) if spacers else '', 'max_spacer_bp':max(spacers) if spacers else '',
                    'spacer_cv':cv,'array_span_bp':starts[-1]+len(motif)-starts[0],
                    'strict_published_art_geometry':strict_art,'broad_art_like_geometry':broad_art,
                })
    # Deduplicate overlapping seed-derived representations of the same array.
    # Prefer more copies, longer motif, lower spacer CV, then lexical motif.
    raw.sort(key=lambda r:(-r['copy_n'],-r['motif_len'],r['spacer_cv'],r['motif']))
    keep=[]
    for r in raw:
        s=set(map(int,r['local_repeat_starts'].split(',')))
        duplicate=False
        for q in keep:
            t=set(map(int,q['local_repeat_starts'].split(',')))
            inter=len(s&t); denom=min(len(s),len(t))
            if denom and inter/denom>=0.8 and abs(r['motif_len']-q['motif_len'])<=12:
                duplicate=True; break
        if not duplicate: keep.append(r)
    return keep

def boolv(x): return str(x).lower()=='true'

def main():
    if len(sys.argv)!=7:
        raise SystemExit('usage: repeat_array_scan.py B4B0_ALL GENOMEDB GENOMEDB_INDEX OUTDIR LOCK_JSON FAMILY_SUMMARY')
    allp,gdb,gidx,outdir,lockp,famsump=sys.argv[1:]
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    lock=json.load(open(lockp)); assert lock['gate_pass'] is True
    rows=list(csv.DictReader(open(allp,encoding='utf-8'),delimiter='\t'))
    if len(rows)!=1605: raise SystemExit(f'FAIL expected 1605 input rows got {len(rows)}')
    if any(r['placement_status']!='UNIQUE_EXACT' for r in rows): raise SystemExit('FAIL non-unique B4B0 row encountered')
    famsum=list(csv.DictReader(open(famsump,encoding='utf-8'),delimiter='\t'))
    fam_order=[r['family_rep_id'] for r in famsum]
    if len(fam_order)!=16: raise SystemExit(f'FAIL family n {len(fam_order)}')

    contigkeys={r['parent_contig_key'] for r in rows}
    gi=load_index(gidx,contigkeys)
    seqs={}
    with open(gdb,'rb',buffering=0) as fh:
        for k,(o,l) in gi.items(): seqs[k]=fetch_db(fh,o,l)

    arrays=[]; locus=[]
    for ix,r in enumerate(rows,1):
        seq=seqs[r['parent_contig_key']]
        s=int(r['candidate_start0']); e=int(r['candidate_end0']); strand=r['strand']
        if strand=='+':
            up=seq[max(0,s-SCAN_FLANK):s]
            dn=seq[e:min(len(seq),e+SCAN_FLANK)]
        elif strand=='-':
            # Candidate-oriented: genomic right side is upstream; genomic left side is downstream.
            up=revcomp(seq[e:min(len(seq),e+SCAN_FLANK)])
            dn=revcomp(seq[max(0,s-SCAN_FLANK):s])
        else: raise SystemExit(f'FAIL strand {strand!r}')
        up_av=len(up); dn_av=len(dn)
        # candidate-relative coords: upstream interval = [-up_av,0), downstream starts at coding length.
        up_arr=scan_region(up,'UPSTREAM',-up_av)
        dn_arr=scan_region(dn,'DOWNSTREAM',int(r['aa_length'])*3)
        arr=up_arr+dn_arr
        for j,a in enumerate(arr,1):
            z={
              'family_rep_id':r['family_rep_id'],'survivor_tier':r['survivor_tier'],'candidate_id':r['candidate_id'],
              'accession':r['accession'],'parent_contig_key':r['parent_contig_key'],'strand':strand,
              'upstream_available_bp':up_av,'downstream_available_bp':dn_av,'array_rank_within_locus':j,
            }; z.update(a); arrays.append(z)
        # best array objective ranking: copies, motif length, regularity, span, upstream before downstream, motif lexical.
        best=None
        if arr:
            best=sorted(arr,key=lambda a:(-a['copy_n'],-a['motif_len'],a['spacer_cv'],-a['array_span_bp'],0 if a['region']=='UPSTREAM' else 1,a['motif']))[0]
        locus.append({
          'family_rep_id':r['family_rep_id'],'survivor_tier':r['survivor_tier'],'candidate_id':r['candidate_id'],'accession':r['accession'],
          'parent_contig_key':r['parent_contig_key'],'strand':strand,'upstream_available_bp':up_av,'downstream_available_bp':dn_av,
          'upstream_complete_500bp':up_av>=500,'upstream_complete_2000bp':up_av>=2000,'upstream_complete_5000bp':up_av>=5000,
          'any_generic_repeat_array':bool(arr),'upstream_repeat_array_n':len(up_arr),'downstream_repeat_array_n':len(dn_arr),
          'strict_published_art_array_n':sum(bool(a['strict_published_art_geometry']) for a in arr),
          'broad_art_like_array_n':sum(bool(a['broad_art_like_geometry']) for a in arr),
          'best_region':best['region'] if best else '','best_motif':best['motif'] if best else '',
          'best_motif_len':best['motif_len'] if best else '','best_copy_n':best['copy_n'] if best else '',
          'best_median_spacer_bp':best['median_spacer_bp'] if best else '',
          'best_strict_published_art_geometry':best['strict_published_art_geometry'] if best else False,
          'best_broad_art_like_geometry':best['broad_art_like_geometry'] if best else False,
          'no_array_in_complete_5kb_upstream':bool(up_av>=ABSENCE_UPSTREAM_COMPLETE and not up_arr),
        })
        if ix%100==0: print(f'[SCAN] {ix}/{len(rows)} arrays={len(arrays)}',flush=True)

    # write arrays and locus table
    arrp=out/'F3B4B1_ALL_DETECTED_REPEAT_ARRAYS.tsv'
    arrcols=list(arrays[0]) if arrays else ['family_rep_id','survivor_tier','candidate_id','accession','parent_contig_key','strand','upstream_available_bp','downstream_available_bp','array_rank_within_locus','region','seed_k','motif','motif_len','copy_n','local_repeat_starts','candidate_relative_repeat_starts','start_to_start_gaps','spacer_lengths','median_spacer_bp','min_spacer_bp','max_spacer_bp','spacer_cv','array_span_bp','strict_published_art_geometry','broad_art_like_geometry']
    with arrp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=arrcols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(arrays)
    locp=out/'F3B4B1_LOCUS_REPEAT_ARRAY_SUMMARY_1605.tsv'; loccols=list(locus[0])
    with locp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=loccols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(locus)

    # family summary: all 1605 scanned, presence may be observed on partial contexts; absence only summarized for complete 5kb upstream instances.
    fs=[]
    for fam in fam_order:
        lr=[x for x in locus if x['family_rep_id']==fam]
        ar=[x for x in arrays if x['family_rep_id']==fam]
        strict=[x for x in ar if bool(x['strict_published_art_geometry'])]
        broad=[x for x in ar if bool(x['broad_art_like_geometry'])]
        up=[x for x in ar if x['region']=='UPSTREAM']
        best=None
        if ar: best=sorted(ar,key=lambda a:(-a['copy_n'],-a['motif_len'],a['spacer_cv'],-a['array_span_bp'],0 if a['region']=='UPSTREAM' else 1,a['motif'],a['candidate_id']))[0]
        complete=[x for x in lr if x['upstream_complete_5000bp']]
        sm={
          'family_rep_id':fam,'survivor_tier':lr[0]['survivor_tier'],'family_instance_n':len(lr),'family_accession_n':len({x['accession'] for x in lr}),
          'instance_n_with_any_repeat_array':sum(x['any_generic_repeat_array'] for x in lr),
          'accession_n_with_any_repeat_array':len({x['accession'] for x in lr if x['any_generic_repeat_array']}),
          'instance_n_with_upstream_repeat_array':len({x['candidate_id'] for x in up}),
          'accession_n_with_upstream_repeat_array':len({x['accession'] for x in up}),
          'strict_published_art_array_n':len(strict),'strict_published_art_instance_n':len({x['candidate_id'] for x in strict}),
          'strict_published_art_accession_n':len({x['accession'] for x in strict}),
          'broad_art_like_array_n':len(broad),'broad_art_like_instance_n':len({x['candidate_id'] for x in broad}),
          'broad_art_like_accession_n':len({x['accession'] for x in broad}),
          'complete_5kb_upstream_instance_n':len(complete),
          'complete_5kb_upstream_no_array_instance_n':sum(x['no_array_in_complete_5kb_upstream'] for x in complete),
          'best_array_candidate_id':best['candidate_id'] if best else '','best_array_accession':best['accession'] if best else '',
          'best_array_region':best['region'] if best else '','best_array_motif':best['motif'] if best else '',
          'best_array_motif_len':best['motif_len'] if best else '','best_array_copy_n':best['copy_n'] if best else '',
          'best_array_median_spacer_bp':best['median_spacer_bp'] if best else '',
          'best_array_strict_published_art_geometry':best['strict_published_art_geometry'] if best else False,
          'best_array_broad_art_like_geometry':best['broad_art_like_geometry'] if best else False,
          'architecture_absence_claim_authorized':False,
        }
        fs.append(sm)
    fsp=out/'F3B4B1_FAMILY_REPEAT_ARRAY_SUMMARY_16.tsv'; fscols=list(fs[0])
    with fsp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fscols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(fs)

    gate={
      'stage':'Stage05F3B4B1','gate':'PASS_F3B4B1_FAMILY_WIDE_REPEAT_ARRAY_ARCHITECTURE_SCAN_COMPLETE','gate_pass':True,
      'family_n':len(fs),'locus_n':len(locus),'detected_repeat_array_n':len(arrays),
      'family_n_with_any_repeat_array':sum(x['instance_n_with_any_repeat_array']>0 for x in fs),
      'family_n_with_upstream_repeat_array':sum(x['instance_n_with_upstream_repeat_array']>0 for x in fs),
      'family_n_with_strict_published_art_geometry':sum(x['strict_published_art_array_n']>0 for x in fs),
      'family_n_with_broad_art_like_geometry':sum(x['broad_art_like_array_n']>0 for x in fs),
      'repeat_detection_role':'prospective generic short-repeat array scan across all 1605 frozen-family loci; ART-paper geometry is descriptive annotation only and not a filter',
      'seed_k_values':list(SEED_K),'min_copies':MIN_COPIES,'start_to_start_gap_range_bp':[GAP_MIN,GAP_MAX],'max_exact_core_motif_bp':MAX_MOTIF,
      'published_art_geometry_annotation':'UPSTREAM; 3-21 copies; every inferred spacer 120-220 bp',
      'broad_art_like_geometry_annotation':'UPSTREAM; 3-30 copies; median spacer 100-240 bp; >=75% spacers 80-280 bp',
      'absence_claim_upstream_completeness_requirement_bp':ABSENCE_UPSTREAM_COMPLETE,
      'architecture_absence_claim_authorized':False,'partner_analysis_performed':False,'novelty_decision_performed':False,
      'all_arrays_sha256':sha256_file(arrp),'locus_summary_sha256':sha256_file(locp),'family_summary_sha256':sha256_file(fsp),
      'claim_boundary':'Repeat-array discovery/geometry only. Presence on partial context is reportable; absence is not a family-level negative claim. Partner gene, ART membership, biochemical RT activity, and novelty remain downstream.'}
    gp=out/'STAGE05F3B4B1_GATE.json'; gp.write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2))
    print('\n=== FAMILY REPEAT ARRAY SUMMARY 16 ===')
    print('\t'.join(fscols))
    for x in fs: print('\t'.join(str(x[c]) for c in fscols))
    print('PASS_STAGE05F3B4B1_FAMILY_WIDE_REPEAT_ARRAY_ARCHITECTURE_SCAN_COMPLETE')

if __name__=='__main__': main()
