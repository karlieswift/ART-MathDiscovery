#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, re, sqlite3, subprocess, sys, tempfile
from pathlib import Path
from collections import defaultdict, Counter

K=31
OV=K-1
LINK_RE=re.compile(r'(?:^|\s)L:([+-]):([^:\s]+):([+-])(?=\s|$)')
RCMAP=str.maketrans('ACGTNacgtn','TGCANtgcan')

def rc(s): return s.translate(RCMAP)[::-1].upper()
def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def read_one_fasta(path):
    name=None; seq=[]
    with open(path) as f:
        for ln in f:
            ln=ln.strip()
            if not ln: continue
            if ln.startswith('>'):
                if name is not None: raise RuntimeError(f'multi-record FASTA: {path}')
                name=ln[1:]
            else: seq.append(ln)
    if name is None: raise RuntimeError(f'empty FASTA: {path}')
    return name,''.join(seq).upper()
def fasta_stream_zst(path, want_seq=True):
    p=subprocess.Popen(['zstdcat',str(path)],stdout=subprocess.PIPE,text=True,bufsize=1)
    hdr=None; seq=[]; seqlen=0
    for line in p.stdout:
        line=line.rstrip('\n\r')
        if not line: continue
        if line.startswith('>'):
            if hdr is not None:
                yield hdr, (''.join(seq).upper() if want_seq else seqlen)
            hdr=line[1:]; seq=[]; seqlen=0
        else:
            s=line.strip()
            if want_seq: seq.append(s)
            else: seqlen += len(s)
    if hdr is not None:
        yield hdr, (''.join(seq).upper() if want_seq else seqlen)
    ret=p.wait()
    if ret: raise RuntimeError(f'zstdcat failed {path} rc={ret}')
def norm_target(acc, tok):
    if tok.startswith(acc+'_'): return tok
    if tok.isdigit(): return f'{acc}_{tok}'
    return tok

def build_meta_db(acc,zst,dbpath):
    if dbpath.exists(): dbpath.unlink()
    cx=sqlite3.connect(dbpath)
    cx.execute('PRAGMA journal_mode=OFF')
    cx.execute('PRAGMA synchronous=OFF')
    cx.execute('PRAGMA temp_store=FILE')
    cx.executescript('CREATE TABLE nodes(id TEXT PRIMARY KEY, len INTEGER NOT NULL); CREATE TABLE edges(src_id TEXT, src_o TEXT, dst_id TEXT, dst_o TEXT);')
    nodes=[]; edges=[]; nnode=nedge=0
    for hdr,seqlen in fasta_stream_zst(zst,want_seq=False):
        uid=hdr.split()[0]
        nodes.append((uid,int(seqlen))); nnode+=1
        for m in LINK_RE.finditer(' '+hdr):
            so,tok,do=m.groups(); dst=norm_target(acc,tok)
            edges.append((uid,so,dst,do)); nedge+=1
        if len(nodes)>=5000:
            cx.executemany('INSERT INTO nodes VALUES (?,?)',nodes); nodes=[]
        if len(edges)>=20000:
            cx.executemany('INSERT INTO edges VALUES (?,?,?,?)',edges); edges=[]
    if nodes: cx.executemany('INSERT INTO nodes VALUES (?,?)',nodes)
    if edges: cx.executemany('INSERT INTO edges VALUES (?,?,?,?)',edges)
    cx.commit()
    cx.execute('CREATE INDEX edge_out ON edges(src_id,src_o)')
    cx.execute('CREATE INDEX edge_in ON edges(dst_id,dst_o)')
    cx.commit()
    return cx,nnode,nedge

def parse_seed(row,acc):
    if row.get('anchor_resolution')!='UNIQUE': raise RuntimeError('non-UNIQUE anchor')
    hdr=row['match_headers'].split(' || ')[0]
    uid=hdr.split()[0]
    orient=row['match_orientations'].split(';')[0]
    pos=int(row['match_positions0'].split(';')[0]); ulen=int(row['match_unitig_lengths'].split(';')[0])
    state_o='+' if orient=='FWD' else '-'
    return uid,state_o,pos,ulen,hdr

def qlen(cx,uid):
    x=cx.execute('SELECT len FROM nodes WHERE id=?',(uid,)).fetchone()
    return None if x is None else int(x[0])
def neighbors(cx,state,direction):
    uid,o=state
    if direction=='RIGHT':
        rows=cx.execute('SELECT DISTINCT dst_id,dst_o FROM edges WHERE src_id=? AND src_o=?',(uid,o)).fetchall()
    else:
        rows=cx.execute('SELECT DISTINCT src_id,src_o FROM edges WHERE dst_id=? AND dst_o=?',(uid,o)).fetchall()
    return [(a,b) for a,b in rows]

def walk(cx,seed_state,direction,need):
    cur=seed_state; visited={cur}; path=[]; recovered=0; stop='ENOUGH'; branch_degree=0
    while recovered < need:
        ns=neighbors(cx,cur,direction); branch_degree=len(ns)
        if len(ns)==0: stop='DEAD_END'; break
        if len(ns)>1: stop='BRANCH'; break
        nxt=ns[0]
        if nxt in visited: stop='CYCLE'; break
        L=qlen(cx,nxt[0])
        if L is None: stop='LINK_TARGET_MISSING'; break
        if L < K: stop='UNITIG_LT_K'; break
        added=L-OV
        path.append((nxt,L,added))
        recovered += added
        visited.add(nxt); cur=nxt
    status='FULL_UNIQUE_PATH' if recovered>=need else 'PARTIAL_GRAPH_CENSORED'
    return {'status':status,'stop_reason':stop,'recovered_bp':recovered,'path':path,'terminal_state':cur,'terminal_degree':branch_degree}

def extract_sequences(zst,wanted):
    out={}
    want=set(wanted)
    for hdr,seq in fasta_stream_zst(zst,want_seq=True):
        uid=hdr.split()[0]
        if uid in want:
            out[uid]=(hdr,seq)
    miss=want-set(out)
    if miss: raise RuntimeError(f'missing wanted unitigs: {sorted(miss)[:10]} n={len(miss)}')
    return out

def oseq(seq,o): return seq if o=='+' else rc(seq)

def build_extension(anchor,walkres,seqs,direction):
    seed_id,seed_o=anchor
    sseq=oseq(seqs[seed_id][1],seed_o)
    ext=''; current_seq=sseq
    # path stored nearest -> farther in direction of traversal
    if direction=='RIGHT':
        prev_state=anchor
        for (state,L,added) in walkres['path']:
            uid,o=state; ns=oseq(seqs[uid][1],o)
            if current_seq[-OV:] != ns[:OV]:
                raise RuntimeError(f'overlap mismatch RIGHT {prev_state}->{state}')
            ext += ns[OV:]
            current_seq=ns; prev_state=state
    else:
        cur_state=anchor
        for (state,L,added) in walkres['path']:
            uid,o=state; ps=oseq(seqs[uid][1],o)
            if ps[-OV:] != current_seq[:OV]:
                raise RuntimeError(f'overlap mismatch LEFT {state}->{cur_state}')
            ext = ps[:-OV] + ext
            current_seq=ps; cur_state=state
    return ext

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--b4a',required=True); ap.add_argument('--b4a2',required=True); ap.add_argument('--out',required=True); ap.add_argument('--work',required=True)
    a=ap.parse_args(); b4a=Path(a.b4a); b4a2=Path(a.b4a2); out=Path(a.out); work=Path(a.work)
    out.mkdir(parents=True,exist_ok=True); work.mkdir(parents=True,exist_ok=True); (out/'contexts_graph_recovered').mkdir(exist_ok=True)
    scanp=b4a2/'outputs/F3B4A2_LOGAN_UNITIG_BOUNDARY_ANCHOR_SCAN.tsv'
    unitigdir=b4a2/'external/logan_unitigs'
    locp=b4a/'outputs/F3B4A_SURVIVOR_NATIVE_LOCUS_MAP.tsv'
    scan=list(csv.DictReader(open(scanp),delimiter='\t')); locs=list(csv.DictReader(open(locp),delimiter='\t'))
    assert len(scan)==31 and all(r['anchor_resolution']=='UNIQUE' for r in scan)
    assert len(locs)==16
    byacc=defaultdict(list)
    for r in scan: byacc[r['source_accession']].append(r)
    boundary=[]; extseq={}
    meta_stats=[]
    for acc in sorted(byacc):
        zst=unitigdir/f'{acc}.unitigs.fa.zst'
        if not zst.exists(): raise FileNotFoundError(zst)
        db=work/f'{acc}.meta.sqlite'
        print(f'[INDEX] {acc} {zst.stat().st_size} bytes',flush=True)
        cx,nnode,nedge=build_meta_db(acc,zst,db)
        print(f'[INDEX_PASS] {acc} nodes={nnode} edges={nedge}',flush=True)
        paths=[]; wanted=set()
        for r in byacc[acc]:
            uid,o,pos,ulen,hdr=parse_seed(r,acc)
            L=qlen(cx,uid)
            if L is None: raise RuntimeError(f'seed not in graph {acc} {uid}')
            if L!=ulen: raise RuntimeError(f'seed length mismatch {uid}: scan={ulen} graph={L}')
            side=r['side']; need=int(r['missing_bp']); direction='LEFT' if side=='LEFT' else 'RIGHT'
            wr=walk(cx,(uid,o),direction,need)
            wanted.add(uid)
            for st,_,_ in wr['path']: wanted.add(st[0])
            paths.append((r,(uid,o),wr))
        seqs=extract_sequences(zst,wanted)
        for r,seed,wr in paths:
            uid,o=seed; anchor=r['anchor31'].upper(); ss=oseq(seqs[uid][1],o)
            if r['side']=='LEFT':
                if ss[:K] != anchor: raise RuntimeError(f'seed anchor mismatch LEFT {r["candidate_id"]}')
            else:
                if ss[-K:] != anchor: raise RuntimeError(f'seed anchor mismatch RIGHT {r["candidate_id"]}')
            ext=build_extension(seed,wr,seqs,'LEFT' if r['side']=='LEFT' else 'RIGHT')
            need=int(r['missing_bp'])
            ext_used=(ext[-need:] if r['side']=='LEFT' else ext[:need]) if need else ''
            full=len(ext)>=need
            key=(r['candidate_id'],r['side'])
            extseq[key]=ext_used
            prow=[]
            for st,L,added in wr['path']: prow.append(f'{st[0]}:{st[1]}:{L}:{added}')
            boundary.append({
                'candidate_id':r['candidate_id'],'source_accession':acc,'side':r['side'],'missing_bp':need,
                'seed_unitig_id':uid,'seed_unitig_orient':o,'seed_unitig_len':len(seqs[uid][1]),
                'graph_status':wr['status'],'stop_reason':wr['stop_reason'],'raw_recovered_bp':wr['recovered_bp'],
                'used_recovered_bp':len(ext_used),'path_node_n':len(wr['path']),'terminal_degree':wr['terminal_degree'],
                'path_states':'|'.join(prow),'extension_sha256':hashlib.sha256(ext_used.encode()).hexdigest() if ext_used else ''
            })
        cx.close(); db.unlink(missing_ok=True)
        meta_stats.append({'accession':acc,'node_n':nnode,'edge_n':nedge,'boundary_n':len(byacc[acc])})
        print(f'[ACCESSION_PASS] {acc} boundaries={len(byacc[acc])}',flush=True)
    bp=out/'F3B4A3_BOUNDARY_GRAPH_EXTENSION_31.tsv'
    with open(bp,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(boundary[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(boundary)
    sp=out/'F3B4A3_GRAPH_METADATA_STATS_12.tsv'
    with open(sp,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(meta_stats[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(meta_stats)
    bby={(r['candidate_id'],r['side']):r for r in boundary}
    ctxrows=[]
    for r in locs:
        cid=r['candidate_id']
        ctx=Path(r['genomic_context_fasta'])
        if not ctx.exists(): ctx=b4a/'outputs/contexts_genomic'/f'{cid}.plusminus10kb.genomic.fna'
        _,base=read_one_fasta(ctx)
        lneed=max(0,10000-int(r['context_left_bp'])); rneed=max(0,10000-int(r['context_right_bp']))
        lext=extseq.get((cid,'LEFT'),''); rext=extseq.get((cid,'RIGHT'),'')
        lfull=(lneed==0 or len(lext)==lneed); rfull=(rneed==0 or len(rext)==rneed)
        recovered=lext+base+rext
        fp=out/'contexts_graph_recovered'/f'{cid}.plusminus10kb.graph_unique.fna'
        fp.write_text(f'>{cid}|GRAPH_UNIQUE_CONTEXT|left={len(lext)}|base={len(base)}|right={len(rext)}\n'+ '\n'.join(recovered[i:i+80] for i in range(0,len(recovered),80))+'\n')
        ctxrows.append({
            'shortlist_rank':r['shortlist_rank'],'family_rep_id':r['family_rep_id'],'candidate_id':cid,'source_accession':r['source_accession'],
            'f3b3d_final_best_tier':r['f3b3d_final_best_tier'],'base_context_bp':len(base),'left_missing_bp':lneed,'right_missing_bp':rneed,
            'left_graph_recovered_bp':len(lext),'right_graph_recovered_bp':len(rext),'left_full_unique':str(lfull),'right_full_unique':str(rfull),
            'full_plusminus10kb_unique_graph_context':str(lfull and rfull),'recovered_context_bp':len(recovered),'recovered_context_fasta':str(fp),
            'recovered_context_sha256':sha256_file(fp),
            'left_stop_reason':('NOT_NEEDED' if lneed==0 else bby[(cid,'LEFT')]['stop_reason']),
            'right_stop_reason':('NOT_NEEDED' if rneed==0 else bby[(cid,'RIGHT')]['stop_reason'])
        })
    cp=out/'F3B4A3_RECOVERED_CONTEXT_16.tsv'
    with open(cp,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(ctxrows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(ctxrows)
    fulln=sum(x['full_plusminus10kb_unique_graph_context']=='True' for x in ctxrows)
    gate={
      'stage':'Stage05F3B4A3','gate':'PASS_F3B4A3_FROZEN_GRAPH_UNIQUE_PATH_CONTEXT_RECOVERY_COMPLETE','gate_pass':True,
      'survivor_n':16,'boundary_anchor_n':31,'unique_anchor_input_n':31,'graph_full_boundary_n':sum(x['graph_status']=='FULL_UNIQUE_PATH' for x in boundary),
      'graph_partial_censored_boundary_n':sum(x['graph_status']!='FULL_UNIQUE_PATH' for x in boundary),'full_plusminus10kb_unique_graph_context_n':fulln,
      'partial_context_n':16-fulln,'kmer_k':K,'unitig_overlap_bp':OV,
      'path_rule':'extend only through degree-1 oriented graph continuation; stop at branch, dead end, cycle, missing link, or after required flank recovered; no abundance-based or annotation-based branch choice',
      'repeat_or_partner_analysis_performed':False,'novelty_decision_performed':False,'architecture_absence_claim_authorized_only_for_full_context_rows':True,
      'boundary_table_sha256':sha256_file(bp),'context_table_sha256':sha256_file(cp),'metadata_stats_sha256':sha256_file(sp),
      'claim_boundary':'Source-context recovery only. Graph branches are never resolved by post-hoc choice. Repeat/partner/ART architecture analysis remains downstream.'
    }
    gp=out/'STAGE05F3B4A3_GATE.json'; gp.write_text(json.dumps(gate,indent=2)+'\n')
    with open(out/'STAGE05F3B4A3_SHA256.tsv','w') as f:
        for p in sorted(out.glob('*')):
            if p.is_file() and p.name!='STAGE05F3B4A3_SHA256.tsv': f.write(f'{sha256_file(p)}\t{p.name}\n')
    print(json.dumps(gate,indent=2))
    print('\n=== GRAPH CONTEXT SUMMARY 16 ===')
    cols=['shortlist_rank','candidate_id','source_accession','f3b3d_final_best_tier','left_graph_recovered_bp','right_graph_recovered_bp','full_plusminus10kb_unique_graph_context','left_stop_reason','right_stop_reason']
    print('\t'.join(cols))
    for x in ctxrows: print('\t'.join(str(x[c]) for c in cols))
    print('PASS_STAGE05F3B4A3_GRAPH_UNIQUE_PATH_CONTEXT_RECOVERY_COMPLETE')
if __name__=='__main__': main()
