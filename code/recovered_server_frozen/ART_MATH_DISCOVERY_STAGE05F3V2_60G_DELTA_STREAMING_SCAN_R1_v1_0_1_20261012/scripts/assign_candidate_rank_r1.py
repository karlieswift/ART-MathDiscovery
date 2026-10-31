#!/usr/bin/env python3
import csv,sys
src,out=sys.argv[1:3]
with open(src,newline='') as fi, open(out,'w',newline='') as fo:
    w=csv.writer(fo,delimiter='\t')
    w.writerow(['stratum','deterministic_rank','accession','selection_hash64','assay_type','organism','tax_id','scientific_name'])
    for rank,line in enumerate(fi,1):
        z=line.rstrip('\n').split('\t')
        if len(z)<7: raise SystemExit(f'FATAL malformed sorted candidate row {rank}: {line[:200]!r}')
        s,h,acc,assay,org,tax,sci=z[:7]
        w.writerow([s,rank,acc,h,assay,org,tax,sci])
