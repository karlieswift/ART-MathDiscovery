#!/usr/bin/env python3
import argparse,csv,hashlib,heapq,re,sys

STRATA=("virome","marine","freshwater","soil","wastewater","human","livestock","other_metagenome")
ALLOWED={"WGS","WGA","Synthetic-Long-Read","WCS","Hi-C","CLONE","CLONEEND","POOLCLONE"}

def norm(s): return (s or '').strip()
def low(s): return norm(s).lower()

def classify(organism, scientific):
    s=(low(organism)+' '+low(scientific)).replace('_',' ')
    if any(x in s for x in ('virome','viral metagenome','virus metagenome','phage metagenome')): return 'virome'
    if any(x in s for x in ('marine metagenome','seawater metagenome','sea water metagenome','ocean metagenome')): return 'marine'
    if any(x in s for x in ('freshwater metagenome','lake water metagenome','groundwater metagenome','river metagenome')): return 'freshwater'
    if any(x in s for x in ('soil metagenome','rhizosphere metagenome')): return 'soil'
    if any(x in s for x in ('wastewater metagenome','activated sludge metagenome','sewage metagenome')): return 'wastewater'
    if any(x in s for x in ('human gut metagenome','human metagenome','human oral metagenome','human skin metagenome','human feces metagenome','human faeces metagenome','human vaginal metagenome','human saliva metagenome','human lung metagenome','human milk metagenome')): return 'human'
    if any(x in s for x in ('bovine metagenome','bovine gut metagenome','pig metagenome','pig gut metagenome','chicken metagenome','chicken gut metagenome','sheep metagenome','sheep gut metagenome','livestock metagenome')): return 'livestock'
    if any(x in s for x in ('metagenome','microbiome','environmental sample')): return 'other_metagenome'
    return None

def pick_col(cols, names):
    by={c.lower():c for c in cols}
    for n in names:
        if n.lower() in by:return by[n.lower()]
    return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',required=True)
    ap.add_argument('--salt',required=True)
    ap.add_argument('--per-stratum',type=int,default=1000)
    a=ap.parse_args()
    r=csv.DictReader(sys.stdin)
    if not r.fieldnames: raise SystemExit('metadata CSV missing header')
    accc=pick_col(r.fieldnames,['acc','accession','run_accession'])
    assayc=pick_col(r.fieldnames,['assay_type','librarystrategy','library_strategy'])
    orgc=pick_col(r.fieldnames,['organism','organism_name'])
    scic=pick_col(r.fieldnames,['scientific_name','sci_name'])
    if not accc or not assayc or not orgc:
        raise SystemExit(f'cannot resolve required columns from {r.fieldnames}')
    heaps={s:[] for s in STRATA}; seen=set(); n=0
    for row in r:
        acc=norm(row.get(accc,'')); assay=norm(row.get(assayc,'')); org=norm(row.get(orgc,'')); sci=norm(row.get(scic,'')) if scic else ''
        if not re.fullmatch(r'[SED]RR\d+',acc): continue
        if assay not in ALLOWED: continue
        if acc in seen: continue
        st=classify(org,sci)
        if st is None: continue
        seen.add(acc); n+=1
        h=hashlib.sha256((a.salt+'|'+acc).encode()).hexdigest(); key=int(h,16)
        item=(-key,acc,assay,org,sci,h)
        q=heaps[st]
        if len(q)<a.per_stratum: heapq.heappush(q,item)
        elif item>q[0]: heapq.heapreplace(q,item)
    rows=[]
    for st in STRATA:
        vals=[(-x[0],)+x[1:] for x in heaps[st]]
        vals.sort(key=lambda x:(x[0],x[1]))
        for rank,x in enumerate(vals,1):
            key,acc,assay,org,sci,h=x
            rows.append({'stratum':st,'stable_rank':rank,'accession':acc,'assay_type':assay,'organism':org,'scientific_name':sci,'stable_sha256':h})
    with open(a.out,'w',newline='') as f:
        cols=['stratum','stable_rank','accession','assay_type','organism','scientific_name','stable_sha256']
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
    sys.stderr.write(f'[PASS] metadata-selected candidate accessions={len(rows)} source-eligible-rows={n}\n')
    for st in STRATA: sys.stderr.write(f'  {st}: {sum(1 for x in rows if x["stratum"]==st)}\n')
if __name__=='__main__':main()
