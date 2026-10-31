#!/usr/bin/env python3
import csv,sys,hashlib,heapq,re,argparse

def norm(x): return (x or '').strip()

def field(row,*names):
    low={k.lower():k for k in row}
    for n in names:
        if n.lower() in low:
            return norm(row[low[n.lower()]])
    return ''

EXTREME=['hydrothermal','hot spring','hotspring','hypersaline','saltern','permafrost','deep sea','deep-sea','hadal','acid mine','acidic','alkaline','geothermal','polar','antarctic','arctic','subseafloor']
HOST=['gut','fecal','faecal','stool','human','oral','saliva','skin','rumen','intestinal','mouse','murine','bovine','cattle','pig','swine','chicken','avian','canine','dog','cat','feline','plant rhizosphere','rhizosphere']
VIRAL=['virome','viral','virus','phage','bacteriophage','metavirome']
ENV=['metagenome','marine','ocean','seawater','freshwater','lake','river','soil','sediment','wastewater','sludge','groundwater','environmental','compost']
BAD_ASSAY=['rna','transcript','chip','atac','bisulfite','methyl','mirna','small rna','est','cloneend']

def stratum(row):
    assay=field(row,'assay_type','assaytype','strategy')
    org=' '.join([field(row,'organism'),field(row,'scientific_name'),field(row,'sci_name')]).lower()
    a=assay.lower()
    if any(x in a for x in BAD_ASSAY): return None
    # Viral first unless clearly extreme/host metagenome; extreme has priority to preserve rare habitat stratum.
    if any(x in org for x in EXTREME): return 'D'
    if any(x in org for x in VIRAL): return 'A'
    is_meta=('meta' in a or 'metagenom' in org)
    if is_meta and any(x in org for x in HOST): return 'C'
    if is_meta or any(x in org for x in ENV): return 'B'
    return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--exclude',required=True)
    ap.add_argument('--k',type=int,default=5000)
    ap.add_argument('--seed',default='F3V2_60G_DELTA_20261011')
    ap.add_argument('--out',required=True)
    a=ap.parse_args()
    exclude={x.strip() for x in open(a.exclude) if x.strip() and not x.startswith('#')}
    heaps={s:[] for s in 'ABCD'}
    rdr=csv.DictReader(sys.stdin)
    if rdr.fieldnames is None: raise SystemExit('metadata CSV has no header')
    n=0; kept=0
    for row in rdr:
        n+=1
        acc=field(row,'acc','accession','run_accession')
        if not acc or acc in exclude: continue
        if not re.fullmatch(r'[SED]RR\d+',acc): continue
        s=stratum(row)
        if s is None: continue
        h=int(hashlib.sha256((a.seed+'|'+s+'|'+acc).encode()).hexdigest()[:16],16)
        # keep K smallest hash values using max-heap encoded as negative
        item=(-h,acc,field(row,'assay_type'),field(row,'organism'),field(row,'tax_id'),field(row,'scientific_name','sci_name'))
        hp=heaps[s]
        if len(hp)<a.k:
            heapq.heappush(hp,item); kept+=1
        elif item[0] > hp[0][0]:
            heapq.heapreplace(hp,item)
    rows=[]
    for s,hp in heaps.items():
        vals=[(-neg,acc,assay,org,tax,sci) for neg,acc,assay,org,tax,sci in hp]
        vals.sort(key=lambda x:(x[0],x[1]))
        for rank,(hv,acc,assay,org,tax,sci) in enumerate(vals,1):
            rows.append((s,rank,acc,hv,assay,org,tax,sci))
    with open(a.out,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['stratum','deterministic_rank','accession','selection_hash64','assay_type','organism','tax_id','scientific_name'])
        w.writerows(rows)
    print(f'metadata_rows={n}',file=sys.stderr)
    for s in 'ABCD': print(f'stratum_{s}_candidate_n={len(heaps[s])}',file=sys.stderr)
    print(f'exclude_n={len(exclude)}',file=sys.stderr)
if __name__=='__main__': main()
