#!/usr/bin/env python3
import csv,hashlib,json,sys,tempfile
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'scripts'))
from common import *
# exact source hashes
with open(ROOT/'config/FROZEN_SOURCE_SHA256.tsv') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        p=ROOT/'frozen_sources'/r['filename']; assert sha256_file(p)==r['sha256'],(r['filename'],sha256_file(p),r['sha256'])
# cohort fixed and untouched by Stage02 development
co=pd.read_csv(ROOT/'config/VALIDATION_COHORT.tsv',sep='\t'); assert co.accession.tolist()==['MW349129','MW349128','MW218148.1']; assert set(co.stage02_development_used)=={'NO'}
# policy exact
p=read_json(ROOT/'config/ENGINE_POLICY.json'); assert p['window_lengths_bp']==[128,256,512,1024,2048,4096]; assert p['proposal_q_threshold']==0.005; assert p['recurrence']['k_values']==[4,6,8]; assert p['recurrence']['min_lag_bp']==8 and p['recurrence']['max_lag_bp']==512; assert p['union_rule']['top_k_cutoff'] is False
# no MarsHill truth constants in execution scripts
for x in (ROOT/'scripts').glob('*.py'):
    s=x.read_text(); assert '213847' not in s and '215071' not in s and 'MW248466' not in s,x
# ACNV 1-6 dimensional contract
m=load_module('tacnv',ROOT/'frozen_sources/ACNV_SOURCE__legacy_acnv_extractor_frozen.py'); seq='ACGT'*80; n,mu,_=m.calculate_values(seq); vals=list(n.values())+list(mu.values())
for k in range(1,7):
    d=m.calculate_all_NV_Kmer_sliding_window(seq,mu,n,k); vals += [float(d.get(t,0.0)) for t in m.generate_kmers(m.standard_aa,k)]
assert len(vals)==5468 and np.isfinite(vals).all(); assert (92,348,1372,5468)==(8+4+16+64,8+4+16+64+256,8+4+16+64+256+1024,8+4+16+64+256+1024+4096)
# generic recurrence sees a synthetic repeated block without any truth-specific lag band
motif='ACGTGTCAGATCAGTA'; syn=(motif+'GATTACA')*50; w=pd.DataFrame([{'start0':0,'end0':len(syn),'window_len':len(syn)}]); z,k,lag,mc=recurrence_scores(syn,w,(4,6,8),8,512,32,2); assert z[0]>0 and 8<=lag[0]<=512
# strict overlap merge: touching intervals do not merge
s=pd.DataFrame({'start0':[0,10,20],'end0':[10,20,31]}); assert merge_overlapping(s)==[[0,10],[10,20],[20,31]]
# GLMY exact parity
B=load_module('tbuilder',ROOT/'frozen_sources/GLMY_DNA_BUILDER__79A2_build_full_glmy5_fresh_features.py'); C=load_module('tcore',ROOT/'frozen_sources/GLMY_CORE__glmy5_core.py'); D=C.load_digraph_class(ROOT/'frozen_sources/GLMY_DIGRAPH__Digraph_Hu.py'); seq='AATGATGCACATCATGGTGGC'; group=('A','GC'); pm=B.get_position_maps(seq,[1,2]); pos,lens=B.group_positions(pm,group); nodes=sorted(pos); edges=[]
for ui,u in enumerate(nodes):
    for vi,v in enumerate(nodes):
        if ui==vi: continue
        d=B.adjusted_distance(u,v,lens.get(u,1),lens.get(v,1))
        if d<=4: edges.append((ui,vi,d))
a=tuple(B.calculate_glmy5(D,len(nodes),tuple(edges))); b=tuple(C.compute_glmy5_from_count(D,len(nodes),tuple(edges)).vector); assert a==b==(2,0,1,12,19),(a,b)
print('[PASS] Stage03A frozen engine source/policy/math selftest')
