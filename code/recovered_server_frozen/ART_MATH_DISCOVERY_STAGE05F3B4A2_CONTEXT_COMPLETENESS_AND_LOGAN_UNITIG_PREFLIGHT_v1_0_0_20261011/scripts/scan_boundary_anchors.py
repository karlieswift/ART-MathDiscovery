#!/usr/bin/env python3
import argparse,csv,subprocess,sys
from pathlib import Path

def rc(s):
    return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def fasta_stream_from_zst(path):
    p=subprocess.Popen(['zstdcat',str(path)],stdout=subprocess.PIPE,text=True,bufsize=1)
    name=None; seq=[]
    for line in p.stdout:
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: yield name,''.join(seq).upper()
            name=line[1:]; seq=[]
        else: seq.append(line)
    if name is not None: yield name,''.join(seq).upper()
    ret=p.wait()
    if ret: raise RuntimeError(f'zstdcat failed {path} rc={ret}')

ap=argparse.ArgumentParser(); ap.add_argument('--anchors',required=True); ap.add_argument('--unitig-dir',required=True); ap.add_argument('--out',required=True)
a=ap.parse_args()
anchors=list(csv.DictReader(open(a.anchors),delimiter='\t'))
by={}
for x in anchors: by.setdefault(x['source_accession'],[]).append(x)
out=[]
for acc,xs in sorted(by.items()):
    p=Path(a.unitig_dir)/f'{acc}.unitigs.fa.zst'
    if not p.exists() or p.stat().st_size==0:
        for x in xs:
            y=dict(x); y.update({'unitig_file_status':'MISSING','anchor_occurrence_n':'','max_immediate_overhang_bp':'','anchor_resolution':'NOT_SCANNED'})
            out.append(y)
        continue
    hits={ (x['candidate_id'],x['side']):[] for x in xs }
    for hdr,seq in fasta_stream_from_zst(p):
        for x in xs:
            k=x['anchor31']; kr=rc(k); key=(x['candidate_id'],x['side'])
            start=0
            while True:
                pos=seq.find(k,start)
                if pos<0: break
                over= pos if x['side']=='LEFT' else len(seq)-(pos+31)
                hits[key].append((hdr,'FWD',pos,len(seq),over)); start=pos+1
            if kr!=k:
                start=0
                while True:
                    pos=seq.find(kr,start)
                    if pos<0: break
                    over= len(seq)-(pos+31) if x['side']=='LEFT' else pos
                    hits[key].append((hdr,'RC',pos,len(seq),over)); start=pos+1
    for x in xs:
        hs=hits[(x['candidate_id'],x['side'])]
        y=dict(x); y['unitig_file_status']='PRESENT'; y['anchor_occurrence_n']=len(hs)
        y['max_immediate_overhang_bp']=max((h[4] for h in hs),default=0)
        y['anchor_resolution']='UNIQUE' if len(hs)==1 else ('ABSENT' if len(hs)==0 else 'AMBIGUOUS')
        y['match_headers']=' || '.join(h[0] for h in hs[:20])
        y['match_orientations']=';'.join(h[1] for h in hs[:20])
        y['match_positions0']=';'.join(str(h[2]) for h in hs[:20])
        y['match_unitig_lengths']=';'.join(str(h[3]) for h in hs[:20])
        out.append(y)
fields=[]
for y in out:
    for k in y:
        if k not in fields: fields.append(k)
with open(a.out,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(out)
print('anchor_rows=',len(out))
print('unique=',sum(x['anchor_resolution']=='UNIQUE' for x in out),'ambiguous=',sum(x['anchor_resolution']=='AMBIGUOUS' for x in out),'absent=',sum(x['anchor_resolution']=='ABSENT' for x in out),'not_scanned=',sum(x['anchor_resolution']=='NOT_SCANNED' for x in out))
