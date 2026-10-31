#!/usr/bin/env python3
"""Stream all source-blind eligible Logan accessions to per-stratum unsorted TSVs.

R1 intentionally removes the fixed top-K acquisition ceiling. Candidate ordering is
still exactly the v1.0.0 SHA256 order: seed|stratum|accession, first 64 bits.
No RT/candidate/model information is read by this program.
"""
import argparse, csv, hashlib, re, sys
from pathlib import Path

EXTREME=['hydrothermal','hot spring','hotspring','hypersaline','saltern','permafrost','deep sea','deep-sea','hadal','acid mine','acidic','alkaline','geothermal','polar','antarctic','arctic','subseafloor']
HOST=['gut','fecal','faecal','stool','human','oral','saliva','skin','rumen','intestinal','mouse','murine','bovine','cattle','pig','swine','chicken','avian','canine','dog','cat','feline','plant rhizosphere','rhizosphere']
VIRAL=['virome','viral','virus','phage','bacteriophage','metavirome']
ENV=['metagenome','marine','ocean','seawater','freshwater','lake','river','soil','sediment','wastewater','sludge','groundwater','environmental','compost']
BAD_ASSAY=['rna','transcript','chip','atac','bisulfite','methyl','mirna','small rna','est','cloneend']
ACC_RE=re.compile(r'[SED]RR\d+')

def norm(x): return (x or '').strip()
def field(row,*names):
    low={k.lower():k for k in row}
    for n in names:
        if n.lower() in low: return norm(row[low[n.lower()]])
    return ''

def stratum(row):
    assay=field(row,'assay_type','assaytype','strategy')
    org=' '.join([field(row,'organism'),field(row,'scientific_name'),field(row,'sci_name')]).lower()
    a=assay.lower()
    if any(x in a for x in BAD_ASSAY): return None
    if any(x in org for x in EXTREME): return 'D'
    if any(x in org for x in VIRAL): return 'A'
    is_meta=('meta' in a or 'metagenom' in org)
    if is_meta and any(x in org for x in HOST): return 'C'
    if is_meta or any(x in org for x in ENV): return 'B'
    return None

def clean(x):
    return norm(x).replace('\t',' ').replace('\r',' ').replace('\n',' ')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--exclude',required=True)
    ap.add_argument('--out-dir',required=True)
    ap.add_argument('--seed',default='F3V2_60G_DELTA_20261011')
    a=ap.parse_args()
    exclude={x.strip() for x in open(a.exclude) if x.strip() and not x.startswith('#')}
    outdir=Path(a.out_dir); outdir.mkdir(parents=True,exist_ok=True)
    fs={s:open(outdir/f'{s}.unsorted.tsv','w',newline='') for s in 'ABCD'}
    counts={s:0 for s in 'ABCD'}; n=0; excluded=0; ineligible=0
    try:
        rdr=csv.DictReader(sys.stdin)
        if rdr.fieldnames is None: raise SystemExit('FATAL: metadata CSV has no header')
        for row in rdr:
            n+=1
            acc=field(row,'acc','accession','run_accession')
            if not acc or not ACC_RE.fullmatch(acc): ineligible+=1; continue
            if acc in exclude: excluded+=1; continue
            s=stratum(row)
            if s is None: ineligible+=1; continue
            h=int(hashlib.sha256((a.seed+'|'+s+'|'+acc).encode()).hexdigest()[:16],16)
            vals=[s,str(h),acc,clean(field(row,'assay_type')),clean(field(row,'organism')),clean(field(row,'tax_id')),clean(field(row,'scientific_name','sci_name'))]
            fs[s].write('\t'.join(vals)+'\n'); counts[s]+=1
    finally:
        for f in fs.values(): f.close()
    print(f'metadata_rows={n}',file=sys.stderr)
    print(f'excluded_old_accession_rows={excluded}',file=sys.stderr)
    print(f'ineligible_rows={ineligible}',file=sys.stderr)
    for s in 'ABCD': print(f'stratum_{s}_eligible_n={counts[s]}',file=sys.stderr)
    print(f'exclude_n={len(exclude)}',file=sys.stderr)
if __name__=='__main__': main()
