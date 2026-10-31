"""Bounded ungapped repeat-core sensitivity, not an RNA-array classifier.
No RT information or mathematical score is passed to this function.
Every emitted copy is checked against an actual sequence template with an
explicit Hamming mismatch limit. A sampled intact 8-mer is needed; indels and
structure-only conservation are outside the assay. Absent hits are not negative.
"""
from collections import defaultdict
from math import floor
from repeat_evidence import entropy

def hamming_at_most(a,b,maximum):
    if len(a)!=len(b):return None
    n=0
    for x,y in zip(a,b):
        n+=x!=y
        if n>maximum:return None
    return n

def greedy_components(positions,L,maxstep):
    groups=[];selected=[]
    for p in sorted(set(positions)):
        if selected and p-selected[-1]>maxstep:
            groups.append(selected);selected=[]
        if not selected or p-selected[-1]>=L:selected.append(p)
    if selected:groups.append(selected)
    return groups

def find_approximate_groups(seq,policy):
    k=int(policy['seed_nt']);idx=defaultdict(list)
    for i in range(len(seq)-k+1):idx[seq[i:i+k]].append(i)
    raw={}
    for L in policy['lengths_nt']:
        offsets=sorted({0,(L-k)//2,L-k});maxmis=floor(L*policy['max_substitution_fraction']+1e-12)
        visited=set()
        for anchor in range(len(seq)-L+1):
            template=seq[anchor:anchor+L]
            if template in visited:continue
            visited.add(template)
            possible=set()
            for off in offsets:
                for p in idx.get(template[off:off+k],()):
                    st=p-off
                    if 0<=st<=len(seq)-L:possible.add(st)
            if len(possible)<policy['min_copies']:continue
            # At least 3 nonoverlapping potential positions must already be
            # connected by the fixed distance rule; split after sequence tests.
            groups=greedy_components(possible,L,policy['max_adjacent_start_distance_bp'])
            if not any(len(g)>=policy['min_copies'] for g in groups):continue
            matches={}
            for st in possible:
                mm=hamming_at_most(template,seq[st:st+L],maxmis)
                if mm is not None:matches[st]=mm
            for ps in greedy_components(matches,L,policy['max_adjacent_start_distance_bp']):
                if len(ps)<policy['min_copies']:continue
                # Same coordinate representation may have more than one
                # template. Choose the one minimizing mismatch sum, then nt.
                key=(tuple(ps),L);mms=[matches[s] for s in ps]
                c={'starts':tuple(ps),'repeat_len':L,'core':template,'mismatches':mms,'seed_count':1,'min_identity_to_template':1-max(mms)/L,'template_entropy_bits':entropy(template)}
                old=raw.get(key)
                if old is None:raw[key]=c
                else:
                    reps=old['seed_count']+1
                    if (sum(mms),template)<(sum(old['mismatches']),old['core']):raw[key]=c
                    raw[key]['seed_count']=reps
    # Collapse contained shifted representations one copy at a time, retaining
    # differing copy counts and non-contained coordinate arrangements.
    buckets=defaultdict(list);kept=[]
    for c in sorted(raw.values(),key=lambda r:(-r['repeat_len'],-len(r['starts']),r['starts'],r['core'])):
        ps=c['starts'];bucket=(len(ps),tuple(p-ps[0] for p in ps))
        parent=next((o for o in buckets[bucket] if all(a<=b and b+c['repeat_len']<=a+o['repeat_len'] for a,b in zip(o['starts'],ps))),None)
        if parent is not None:parent['seed_count']+=c['seed_count']
        else:buckets[bucket].append(c);kept.append(c)
    for r in kept:
        assert all(hamming_at_most(r['core'],seq[p:p+r['repeat_len']],floor(r['repeat_len']*policy['max_substitution_fraction']+1e-12)) is not None for p in r['starts'])
    return sorted(kept,key=lambda r:(r['starts'][0],r['starts'][-1]+r['repeat_len'],r['core']))
