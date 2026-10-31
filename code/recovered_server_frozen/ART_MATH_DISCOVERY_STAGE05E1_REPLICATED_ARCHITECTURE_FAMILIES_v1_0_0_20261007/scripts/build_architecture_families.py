#!/usr/bin/env python3
import argparse, json, hashlib
from pathlib import Path
from collections import defaultdict, Counter, deque
import pandas as pd

TRANS=str.maketrans('ACGTN','TGCAN')
def canon_nt(s):
    s=str(s).upper().strip()
    rc=s.translate(TRANS)[::-1]
    return min(s,rc)

def semis(vals):
    return ';'.join(sorted({str(x) for x in vals if str(x) not in ('','nan','None')}))

def parse_keys(x):
    if pd.isna(x): return []
    return [z for z in str(x).split(';') if z]

def build_graph(thr, keys, seqmap, cand_core, aa_to_fams, lk):
    cand_fams=defaultdict(set)
    sub=lk[lk['distance_to_math_candidate_bp']<=thr]
    for r in sub.itertuples(index=False):
        for fc in aa_to_fams.get(str(r.aa_sha256),()):
            cand_fams[str(r.candidate_key)].add(fc)
    inv=defaultdict(list)
    for ck,fs in cand_fams.items():
        for fc in fs: inv[fc].append(ck)
    adj=defaultdict(set); ev={}
    for fc, ks0 in inv.items():
        ks=sorted(set(ks0))
        for i,a in enumerate(ks):
            ca=cand_core.get(a,set())
            if not ca: continue
            for b in ks[i+1:]:
                cores=ca & cand_core.get(b,set())
                if cores:
                    adj[a].add(b); adj[b].add(a)
                    k=(a,b)
                    e=ev.setdefault(k,{'fams':set(),'cores':set()})
                    e['fams'].add(fc); e['cores'].update(cores)
    comps=[]; seen=set()
    for k in sorted(adj):
        if k in seen: continue
        q=[k]; seen.add(k); c=[]
        while q:
            u=q.pop(); c.append(u)
            for v in sorted(adj[u]):
                if v not in seen:
                    seen.add(v); q.append(v)
        comps.append(sorted(c))
    repl=[c for c in comps if len({seqmap[x] for x in c})>=2]
    repl=sorted(repl,key=lambda c:(-len({seqmap[x] for x in c}),-len(c),c[0]))
    cid={}
    prefix=f'ARCH{thr}'
    for i,c in enumerate(repl,1):
        fid=f'{prefix}_{i:04d}'
        for x in c: cid[x]=fid
    return ev,repl,cid,cand_fams

def family_table(repl,cid,seqmap,stratmap,rankmap,ev,thr):
    rows=[]
    for comp in repl:
        fid=cid[comp[0]]
        nodes=set(comp)
        ies=[(a,b,e) for (a,b),e in ev.items() if a in nodes and b in nodes]
        fams=set(); cores=set()
        for _,_,e in ies:
            fams |= e['fams']; cores |= e['cores']
        n1=sum(stratmap[x].startswith('N1_') for x in comp)
        n2=sum(stratmap[x].startswith('N2_') for x in comp)
        rows.append({
            'family_id':fid,
            'threshold_bp':thr,
            'candidate_n':len(comp),
            'genome_n':len({seqmap[x] for x in comp}),
            'N1_candidate_n':n1,
            'N2_candidate_n':n2,
            'min_region_rank':min(rankmap[x] for x in comp),
            'edge_n':len(ies),
            'shared_FAM50_cluster_n':len(fams),
            'shared_exact_repeat_core_n':len(cores),
            'candidate_keys':';'.join(comp),
            'sequence_ids':';'.join(sorted({seqmap[x] for x in comp})),
            'shared_FAM50_clusters':';'.join(sorted(fams)),
            'shared_repeat_cores_canonical':';'.join(sorted(cores)),
            'family_review_stratum':'E1A_PRIMARY_REPLICATED_WITH_N1' if thr==250 and n1>0 else ('E1B_PRIMARY_REPLICATED_N2_ONLY' if thr==250 else 'SENSITIVITY_REPLICATED')
        })
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--d1-root',required=True)
    ap.add_argument('--fam50-tsv',required=True)
    ap.add_argument('--outdir',required=True)
    a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    d=Path(a.d1_root)
    can=pd.read_csv(d/'config/NOVELTY_FIRST_999_CANDIDATES.tsv',sep='\t')
    links=pd.read_csv(d/'config/NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv',sep='\t')
    orf=pd.read_csv(d/'config/NOVELTY_FIRST_24020_ORF_REGISTRY.tsv',sep='\t')
    rep=pd.read_csv(d/'config/NOVELTY_FIRST_NONCODING_REPEAT_GROUPS.tsv',sep='\t')
    fam=pd.read_csv(a.fam50_tsv,sep='\t')
    n12=can[can['novelty_architecture_stratum'].astype(str).str.startswith(('N1_','N2_'))].copy()
    if len(n12)!=602: raise SystemExit(f'Expected 602 N1/N2 candidates, got {len(n12)}')
    keys=set(n12.candidate_key.astype(str)); seqmap=n12.set_index('candidate_key').sequence_id.astype(str).to_dict(); stratmap=n12.set_index('candidate_key').novelty_architecture_stratum.astype(str).to_dict(); rankmap=n12.set_index('candidate_key').region_rank.astype(int).to_dict()
    aa_to_fams=defaultdict(set)
    for r in fam.itertuples(index=False): aa_to_fams[str(r.member_aa_sha256)].add(str(r.cluster_id))
    lk=links.merge(orf[['record_key','aa_sha256']],on='record_key',how='left')
    lk=lk[lk.candidate_key.astype(str).isin(keys)].copy(); lk['candidate_key']=lk.candidate_key.astype(str)
    rr=rep[(rep.candidate_key.astype(str).isin(keys))&(rep.high_complexity==1)&(rep.full_span_outside_predicted_CDS==1)&(rep.overlap_candidate_bp>0)].copy()
    rr['candidate_key']=rr.candidate_key.astype(str); rr['core_canonical']=rr.core_nt.map(canon_nt)
    cand_core=rr.groupby('candidate_key')['core_canonical'].apply(set).to_dict()
    ev250, comps250, id250, cf250=build_graph(250,keys,seqmap,cand_core,aa_to_fams,lk)
    ev1000, comps1000, id1000, cf1000=build_graph(1000,keys,seqmap,cand_core,aa_to_fams,lk)
    def edge_df(ev,thr):
        rows=[]
        for (a1,b1),e in sorted(ev.items()):
            rows.append({'threshold_bp':thr,'candidate_a':a1,'candidate_b':b1,'sequence_a':seqmap[a1],'sequence_b':seqmap[b1],
                         'shared_FAM50_cluster_n':len(e['fams']),'shared_exact_repeat_core_n':len(e['cores']),
                         'shared_FAM50_clusters':';'.join(sorted(e['fams'])),'shared_repeat_cores_canonical':';'.join(sorted(e['cores']))})
        return pd.DataFrame(rows)
    edge_df(ev250,250).to_csv(out/'STAGE05E1_PRIMARY_EDGES.tsv',sep='\t',index=False)
    edge_df(ev1000,1000).to_csv(out/'STAGE05E1_SENSITIVITY_EDGES.tsv',sep='\t',index=False)
    ft250=family_table(comps250,id250,seqmap,stratmap,rankmap,ev250,250)
    ft1000=family_table(comps1000,id1000,seqmap,stratmap,rankmap,ev1000,1000)
    ft250.to_csv(out/'STAGE05E1_PRIMARY_REPLICATED_FAMILIES.tsv',sep='\t',index=False)
    ft1000.to_csv(out/'STAGE05E1_SENSITIVITY_REPLICATED_FAMILIES.tsv',sep='\t',index=False)
    rows=[]; counts=Counter()
    for r in n12.sort_values('candidate_key').itertuples(index=False):
        ck=str(r.candidate_key)
        if ck in id250:
            famrow=ft250[ft250.family_id==id250[ck]].iloc[0]
            st='E1A_PRIMARY_REPLICATED_WITH_N1' if int(famrow.N1_candidate_n)>0 else 'E1B_PRIMARY_REPLICATED_N2_ONLY'
        elif ck in id1000: st='E1C_SENSITIVITY_ONLY_REPLICATED'
        else: st='E1D_UNREPLICATED_RETAINED'
        counts[st]+=1
        rows.append({'candidate_key':ck,'sequence_id':str(r.sequence_id),'region_rank':int(r.region_rank),'novelty_architecture_stratum':str(r.novelty_architecture_stratum),
                     'primary_family_id':id250.get(ck,''),'sensitivity_family_id':id1000.get(ck,''),'stage05e1_stratum':st,
                     'qualifying_repeat_core_n':len(cand_core.get(ck,set())),'compact250_FAM50_n':len(cf250.get(ck,set())),'compact1000_FAM50_n':len(cf1000.get(ck,set()))})
    pd.DataFrame(rows).to_csv(out/'STAGE05E1_CANDIDATE_FAMILY_ASSIGNMENTS.tsv',sep='\t',index=False)
    summary={
      'stage':'Stage05E1','n_N1_N2_candidates':602,
      'primary_orf_distance_bp':250,'sensitivity_orf_distance_bp':1000,
      'primary_edge_n':len(ev250),'primary_replicated_family_n':len(comps250),'primary_replicated_candidate_n':len(id250),
      'sensitivity_edge_n':len(ev1000),'sensitivity_replicated_family_n':len(comps1000),'sensitivity_replicated_candidate_n':len(id1000),
      'candidate_stratum_counts':dict(counts),
      'family_construction_uses_structure_results':False,'candidate_filtering_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,
      'gate':'PASS_REPLICATED_ARCHITECTURE_FAMILIES_FROZEN','gate_pass':True,
      'claim_boundary':'Architecture families organize frozen N1/N2 candidate loci by shared repeat-core plus nearby protein-family evidence; they do not establish RT activity or ART identity.'
    }
    (out/'STAGE05E1_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'GATE_STAGE05E1.json').write_text(json.dumps(summary,indent=2)+'\n')
    # sha manifest for outputs except itself
    files=sorted([p for p in out.iterdir() if p.is_file() and p.name!='STAGE05E1_SHA256.tsv'])
    with open(out/'STAGE05E1_SHA256.tsv','w') as h:
        for p in files:
            h.write(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
