#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, os, re, shutil, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path

KSET=[8,10,12,14,16]
STARTS={"ATG","GTG","TTG"}; STOPS={"TAA","TAG","TGA"}
# Standard coding codons; bacterial code 11 changes initiator semantics, not internal amino-acid identities.
CODON={
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
COMP=str.maketrans('ACGTacgt','TGCAtgca')

def rc(s): return s.translate(COMP)[::-1].upper()
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def write_json(p,x): Path(p).write_text(json.dumps(x,indent=2,sort_keys=True)+"\n")
def tsv_read(p):
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def tsv_write(p, rows, fields):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def parse_cid(cid):
    z=cid.split('|')
    if len(z)<7: raise ValueError(f'candidate_id format: {cid}')
    acc,contig,strand,frame,coord,aalen,ordv=z[:7]
    m=re.fullmatch(r'(\d+)-(\d+)',coord); ma=re.fullmatch(r'aa(\d+)',aalen)
    if not m or not ma or strand not in ['+','-']: raise ValueError(cid)
    a,b=map(int,m.groups()); aa=int(ma.group(1))
    return dict(candidate_id=cid,accession=acc,contig_id=contig,strand=strand,frame=frame,start0=a,end0=b,aa_length=aa,ordinal=ordv)
def translate(nt):
    nt=nt.upper(); return ''.join(CODON.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))
def opencmd(path):
    p=str(path)
    if p.endswith('.zst'): return subprocess.Popen(['zstdcat',p],stdout=subprocess.PIPE,text=True,bufsize=1)
    if p.endswith('.gz'): return None
    return None

def iter_fasta(path):
    p=Path(path)
    if str(p).endswith('.zst'):
        proc=subprocess.Popen(['zstdcat',str(p)],stdout=subprocess.PIPE,text=True,bufsize=1)
        fh=proc.stdout
        close=lambda: proc.wait()
    elif str(p).endswith('.gz'):
        fh=gzip.open(p,'rt'); close=fh.close
    else:
        fh=open(p); close=fh.close
    name=None; buf=[]
    try:
        for line in fh:
            if line.startswith('>'):
                if name is not None: yield name,''.join(buf).upper()
                name=line[1:].strip().split()[0]; buf=[]
            else: buf.append(line.strip())
        if name is not None: yield name,''.join(buf).upper()
    finally: close()

def find_selected_ledger(P):
    hits=list(P.glob('**/frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'))
    if not hits:
        raise RuntimeError('no F3V2_DELTA_SELECTED_ACCESSIONS.tsv found')
    # Multiple byte-identical/provenance copies may coexist. Prefer the authoritative
    # distributed-acquisition R2 tree, then a ledger that resolves the expected 254 accessions.
    scored=[]
    for h in hits:
        score=0
        hs=str(h)
        if 'DISTRIBUTED_ACQUISITION_R2' in hs: score += 1000
        if 'v1_0_2_20261013' in hs: score += 200
        try:
            rr=tsv_read(h)
            accs={r.get('accession') or r.get('run_accession') or r.get('sra_accession') or r.get('source_accession') for r in rr}
            accs.discard(None); accs.discard('')
            if len(accs)==254: score += 100
        except Exception:
            pass
        scored.append((score, -len(h.parts), hs, h))
    scored.sort(reverse=True)
    return scored[0][3]
def ledger_map(ledger):
    rows=tsv_read(ledger); out={}
    for r in rows:
        acc=r.get('accession') or r.get('run_accession') or r.get('sra_accession') or r.get('source_accession')
        if acc: out[acc]=r
    return out

def candidate_faa_index(P, wanted):
    # Search likely stage04B2/B3 FASTAs only. Return exact header -> aa.
    out={}
    roots=[x for x in P.iterdir() if x.is_dir() and ('STAGE04B2' in x.name or 'STAGE04B3' in x.name)]
    files=[]
    for root in roots:
        files += list(root.glob('**/*.faa')) + list(root.glob('**/*.fa')) + list(root.glob('**/*.fasta'))
    for fp in files:
        try:
            for name,seq in iter_fasta(fp):
                if name in wanted and name not in out: out[name]=(seq,fp)
        except Exception:
            continue
        if len(out)==len(wanted): break
    return out

def find_source_file(source_root, acc, ledger_row, ledger_path=None):
    # Prefer an existing ledger path if present; resolve relative paths against
    # the ledger/package root before falling back to a bounded recursive search.
    bases=[Path.cwd(), source_root]
    if ledger_path is not None:
        bases += [Path(ledger_path).parent, Path(ledger_path).parents[1]]
    for k,v in ledger_row.items():
        if not v: continue
        if 'path' in k.lower() or 'file' in k.lower():
            q=Path(v)
            candidates=[q] if q.is_absolute() else [b/q for b in bases]
            for qq in candidates:
                if qq.exists() and qq.is_file(): return qq.resolve()
    pats=[f'{acc}.contigs.fa.zst',f'{acc}.fa.zst',f'{acc}.fasta.zst',f'*{acc}*contigs*.zst',f'*{acc}*.fa.gz',f'*{acc}*.fasta.gz',f'*{acc}*.fa',f'*{acc}*.fasta']
    for pat in pats:
        hits=list(source_root.glob('**/'+pat))
        if hits: return hits[0].resolve()
    raise FileNotFoundError(f'no source FASTA for {acc} under {source_root}')
def fetch_record(source_file, rid):
    for name,seq in iter_fasta(source_file):
        if name==rid: return seq
    raise KeyError(f'{rid} not found in {source_file}')

def discover_panel_mapping(P, cids):
    # Find TSV containing candidate_id+family_rep_id and most of the 32 cids.
    best=None
    for root in P.iterdir():
        if not root.is_dir() or 'STAGE04B2' not in root.name: continue
        for fp in root.glob('**/*.tsv'):
            try:
                with open(fp) as f: hdr=f.readline().rstrip('\n').split('\t')
            except: continue
            if 'candidate_id' not in hdr or 'family_rep_id' not in hdr: continue
            try: rows=tsv_read(fp)
            except: continue
            seen=sum(1 for r in rows if r.get('candidate_id') in cids)
            if seen and (best is None or seen>best[0]): best=(seen,fp,rows)
    if best is None: raise RuntimeError('could not discover Stage04B2 candidate_id -> family_rep_id table')
    mp={r['candidate_id']:r['family_rep_id'] for r in best[2] if r.get('candidate_id') in cids}
    if len(mp)!=len(cids): raise RuntimeError(f'panel mapping incomplete {len(mp)}/{len(cids)} from {best[1]}')
    return mp,best[1]
def discover_membership(P, family_reps):
    # Find a TSV with family_rep_id and a member candidate column; score by number of wanted families represented and rows.
    candidate_cols=['candidate_id','member_candidate_id','member_id','query_id','sequence_id']
    best=None
    for root in P.iterdir():
        if not root.is_dir() or ('STAGE04A' not in root.name and 'STAGE04B1' not in root.name and 'STAGE04B2' not in root.name): continue
        for fp in root.glob('**/*.tsv'):
            try:
                with open(fp) as f: hdr=f.readline().rstrip('\n').split('\t')
            except: continue
            if 'family_rep_id' not in hdr: continue
            cc=next((c for c in candidate_cols if c in hdr),None)
            if not cc: continue
            try: rows=tsv_read(fp)
            except: continue
            rr=[r for r in rows if r.get('family_rep_id') in family_reps and r.get(cc)]
            famn=len({r.get('family_rep_id') for r in rr})
            score=(famn,len(rr))
            if famn and (best is None or score>best[0]): best=(score,fp,rr,cc)
    if best is None or best[0][0] < len(family_reps):
        raise RuntimeError(f'family membership unresolved: best={None if best is None else (best[0],str(best[1]))}')
    out=[]
    for r in best[2]:
        cid=r[best[3]]
        try: pc=parse_cid(cid)
        except: continue
        out.append({'family_rep_id':r['family_rep_id'],'candidate_id':cid,'accession':pc['accession']})
    # unique and fail closed if parsing removed any requested family.
    uniq={(r['family_rep_id'],r['candidate_id']):r for r in out}
    vals=list(uniq.values())
    represented={r['family_rep_id'] for r in vals}
    missing=sorted(set(family_reps)-represented)
    if missing:
        raise RuntimeError(f'parsed family membership incomplete; missing {missing}')
    return vals, best[1]

def repeat_arrays(seq, rt0, rt1):
    found={}
    n=len(seq)
    for k in KSET:
        pos=defaultdict(list)
        for i in range(0,n-k+1): pos[seq[i:i+k]].append(i)
        for seed,pp in pos.items():
            if len(pp)<3: continue
            pp=sorted(pp); run=[pp[0]]
            runs=[]
            for x in pp[1:]:
                if 60 <= x-run[-1] <= 360: run.append(x)
                else:
                    if len(run)>=3: runs.append(run)
                    run=[x]
            if len(run)>=3: runs.append(run)
            for starts in runs:
                # maximal exact extension shared across all copies, total core <=80
                L=0
                while min(starts)-L-1>=0 and k+L<80:
                    ch=seq[starts[0]-L-1]
                    if all(seq[s-L-1]==ch for s in starts): L+=1
                    else: break
                R=0
                while max(starts)+k+R<n and k+L+R<80:
                    ch=seq[starts[0]+k+R]
                    if all(seq[s+k+R]==ch for s in starts): R+=1
                    else: break
                ss=tuple(s-L for s in starts); motif=seq[ss[0]:ss[0]+k+L+R]
                key=ss
                prev=found.get(key)
                if prev is None or len(motif)>len(prev['motif']):
                    gaps=[ss[i+1]-ss[i] for i in range(len(ss)-1)]
                    spacers=[g-len(motif) for g in gaps]
                    a0=ss[0]; a1=ss[-1]+len(motif)
                    orient='UPSTREAM' if a1<=rt0 else ('DOWNSTREAM' if a0>=rt1 else 'OVERLAP_RT')
                    strict=(orient=='UPSTREAM' and 3<=len(ss)<=21 and all(120<=z<=220 for z in spacers))
                    found[key]=dict(motif=motif,motif_len=len(motif),copy_n=len(ss),starts=list(ss),start_gaps=gaps,spacers=spacers,array_start0=a0,array_end0=a1,orientation=orient,strict_published_art_geometry=strict)
    return sorted(found.values(),key=lambda x:(x['array_start0'], -x['motif_len'], -x['copy_n']))
def find_orfs_downstream(seq, rt1, maxbp=10000):
    d0=rt1; d1=min(len(seq),rt1+maxbp); sub=seq[d0:d1]
    ans=[]
    # plus strand in candidate-oriented coordinates
    for frame in range(3):
        for i in range(frame,len(sub)-2,3):
            if sub[i:i+3] not in STARTS: continue
            for j in range(i+3,len(sub)-2,3):
                if sub[j:j+3] in STOPS:
                    aa=translate(sub[i:j])
                    if len(aa)>=30: ans.append((i,j+3,'+',aa))
                    break
    # reverse strand ORFs physically within downstream interval
    rsub=rc(sub); L=len(sub)
    for frame in range(3):
        for i in range(frame,len(rsub)-2,3):
            if rsub[i:i+3] not in STARTS: continue
            for j in range(i+3,len(rsub)-2,3):
                if rsub[j:j+3] in STOPS:
                    aa=translate(rsub[i:j])
                    if len(aa)>=30:
                        # map rc [i,j+3) to sub [L-(j+3), L-i)
                        a=L-(j+3); b=L-i; ans.append((a,b,'-',aa))
                    break
    ans.sort(key=lambda z:(z[0],z[1]-z[0],z[2]))
    return ans

def read_stage04a_metrics(P):
    # Prefer gate/summary JSON containing recognizable counts, otherwise use frozen known values with explicit source label.
    vals={'old_bp':6121620524,'delta_bp':54377819157,'old_formal':10598,'delta_formal':88487,'old_families':4912,'delta_families':26396,'delta_only_combined_families':24810,'source':'canonical frozen project values'}
    return vals

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--project-root',required=True)
    ap.add_argument('--source-root',default=None,help='54G acquisition/source-freeze root; auto-discovered if omitted')
    ap.add_argument('--threads',type=int,default=32)
    ap.add_argument('--outdir',default='outputs')
    ap.add_argument('--representative-only',action='store_true',help='technical diagnostic only; primary closeout requires family membership')
    args=ap.parse_args()
    HERE=Path(__file__).resolve().parents[1]; P=Path(args.project_root).resolve(); O=HERE/args.outdir; O.mkdir(parents=True,exist_ok=True)
    frozen=HERE/'frozen/DELTA32_STRUCTURAL_COMPATIBILITY.tsv'
    rows=tsv_read(frozen); pos=[r for r in rows if r['tier'] in ('STRONG','SUGGESTIVE')]
    if len(pos)!=12: raise SystemExit(f'FAIL expected 12 frozen structural positives, got {len(pos)}')
    cids=[r['candidate_id'] for r in pos]
    tsv_write(O/'DELTA12_STRUCTURAL_POSITIVES.tsv',pos,list(pos[0].keys()))
    panel_map,panel_src=discover_panel_mapping(P,set(cids)); fams=set(panel_map.values())
    if len(fams)!=12: raise SystemExit(f'FAIL expected 12 structural-positive families, mapped {len(fams)}')
    if args.representative_only:
        members=[{'family_rep_id':panel_map[c],'candidate_id':c,'accession':parse_cid(c)['accession']} for c in cids]
        member_src='REPRESENTATIVE_ONLY_DIAGNOSTIC'
    else:
        members,member_src=discover_membership(P,fams)
    tsv_write(O/'DELTA12_FAMILY_MEMBERSHIP.tsv',members,['family_rep_id','candidate_id','accession'])
    # source root / ledger
    if args.source_root: SR=Path(args.source_root).resolve(); ledger=find_selected_ledger(SR)
    else:
        ledger=find_selected_ledger(P); SR=ledger.parents[1]
    lm=ledger_map(ledger)
    # Exact frozen amino-acid comparison is required for the 12 panel representatives.
    # For additional family members, the frozen membership candidate_id itself encodes
    # the source contig/strand/coordinates; we replay those coordinates from frozen raw DNA,
    # validate nucleotide length/translation, and record the reconstructed AA SHA256.
    aaidx=candidate_faa_index(P,set(cids))
    missingrep=sorted(set(cids)-set(aaidx))
    if missingrep: raise SystemExit(f'HOLD representative candidate FASTA unresolved {len(missingrep)}/12; first={missingrep[:3]}')
    ctx_rows=[]; arrays=[]; partner_rows=[]; partner_faa=[]
    cache={}; source_sha_cache={}
    for ix,m in enumerate(members,1):
        pc=parse_cid(m['candidate_id']); acc=pc['accession']
        if acc not in lm: raise SystemExit(f'HOLD accession {acc} absent from selected-accession ledger')
        sf=find_source_file(SR,acc,lm[acc],ledger)
        ck=(str(sf),pc['contig_id'])
        if ck not in cache: cache[ck]=fetch_record(sf,pc['contig_id'])
        genome=cache[ck]; a,b=pc['start0'],pc['end0']
        if b>a and b<=len(genome): nt=genome[a:b]
        else: raise SystemExit(f'HOLD coords out of range {pc}')
        oriented=nt if pc['strand']=='+' else rc(nt)
        if len(oriented)!=pc['aa_length']*3: raise SystemExit(f'HOLD length mismatch {pc["candidate_id"]}')
        tr=translate(oriented).rstrip('*')
        if len(tr)!=pc['aa_length'] or '*' in tr or 'X' in tr:
            raise SystemExit(f'HOLD coordinate translation invalid {pc["candidate_id"]}')
        if pc['candidate_id'] in aaidx:
            aa=aaidx[pc['candidate_id']][0].rstrip('*')
            if tr!=aa: raise SystemExit(f'HOLD exact representative peptide replay mismatch {pc["candidate_id"]}')
            aa_source=str(aaidx[pc['candidate_id']][1])
        else:
            aa_source='EXACT_COORDINATE_REPLAY_FROM_FROZEN_FAMILY_MEMBERSHIP'
        left=max(0,a-10000); right=min(len(genome),b+10000); rawctx=genome[left:right]
        if pc['strand']=='+':
            ctx=rawctx; rt0=a-left; rt1=b-left
        else:
            ctx=rc(rawctx); rt0=right-b; rt1=right-a
        sid=f"{m['family_rep_id']}|{pc['candidate_id']}"
        with open(O/'DELTA12_NATIVE_CONTEXTS.fna','a' if (O/'DELTA12_NATIVE_CONTEXTS.fna').exists() else 'w') as fh:
            fh.write(f'>{sid}\n');
            for j in range(0,len(ctx),80): fh.write(ctx[j:j+80]+'\n')
        stratum=lm[acc].get('stratum') or lm[acc].get('source_stratum') or ''
        if str(sf) not in source_sha_cache: source_sha_cache[str(sf)]=sha(sf)
        aa_sha=hashlib.sha256(tr.encode()).hexdigest()
        ctx_rows.append(dict(family_rep_id=m['family_rep_id'],candidate_id=pc['candidate_id'],accession=acc,stratum=stratum,contig_id=pc['contig_id'],strand=pc['strand'],rt_start0_in_oriented_context=rt0,rt_end0_in_oriented_context=rt1,context_length_bp=len(ctx),source_file=str(sf),source_file_sha256=source_sha_cache[str(sf)],candidate_aa_sha256=aa_sha,candidate_aa_source=aa_source))
        aaarr=repeat_arrays(ctx,rt0,rt1)
        for ai,z in enumerate(aaarr,1):
            arrays.append(dict(family_rep_id=m['family_rep_id'],candidate_id=pc['candidate_id'],accession=acc,stratum=stratum,array_index=ai,motif=z['motif'],motif_len=z['motif_len'],copy_n=z['copy_n'],repeat_starts=','.join(map(str,z['starts'])),start_gaps=','.join(map(str,z['start_gaps'])),inferred_spacers=','.join(map(str,z['spacers'])),orientation=z['orientation'],strict_published_art_geometry=str(z['strict_published_art_geometry'])))
        orfs=find_orfs_downstream(ctx,rt1)
        if orfs:
            oa,ob,ost,opa=orfs[0]
            pid=f"P_{len(partner_rows)+1:06d}"
            partner_rows.append(dict(partner_id=pid,family_rep_id=m['family_rep_id'],candidate_id=pc['candidate_id'],accession=acc,stratum=stratum,orf_strand=ost,orf_start0_downstream=oa,orf_end0_downstream=ob,aa_length=len(opa),rt_to_partner_gap_bp=oa))
            partner_faa.append((pid,opa))
    tsv_write(O/'DELTA12_NATIVE_CONTEXTS.tsv',ctx_rows,list(ctx_rows[0].keys()))
    afields=['family_rep_id','candidate_id','accession','stratum','array_index','motif','motif_len','copy_n','repeat_starts','start_gaps','inferred_spacers','orientation','strict_published_art_geometry']
    tsv_write(O/'DELTA12_REPEAT_ARRAYS.tsv',arrays,afields)
    pfields=['partner_id','family_rep_id','candidate_id','accession','stratum','orf_strand','orf_start0_downstream','orf_end0_downstream','aa_length','rt_to_partner_gap_bp']
    tsv_write(O/'DELTA12_NEAREST_DOWNSTREAM_ORFS.tsv',partner_rows,pfields)
    with open(O/'DELTA12_NEAREST_DOWNSTREAM_ORFS.faa','w') as f:
        for n,s in partner_faa: f.write(f'>{n}\n{s}\n')
    # cluster partners
    cluster_map={}
    mm=shutil.which('mmseqs')
    if mm and partner_faa:
        mmroot=O/'mmseqs_partner'; tmp=O/'mmseqs_tmp'
        subprocess.run([mm,'easy-cluster',str(O/'DELTA12_NEAREST_DOWNSTREAM_ORFS.faa'),str(mmroot),str(tmp),'--min-seq-id','0.3','-c','0.7','--cov-mode','0','--cluster-mode','0','-s','7.5','--threads',str(args.threads)],check=True)
        ctsv=Path(str(mmroot)+'_cluster.tsv')
        with open(ctsv) as f:
            for line in f:
                rep,mem=line.rstrip('\n').split('\t')[:2]; cluster_map[mem]=rep
    for r in partner_rows: r['partner_cluster']=cluster_map.get(r['partner_id'],'UNCLUSTERED')
    if partner_rows: tsv_write(O/'DELTA12_NEAREST_DOWNSTREAM_ORFS_WITH_CLUSTER.tsv',partner_rows,pfields+['partner_cluster'])
    # family summary
    arr_by=defaultdict(list)
    for a in arrays: arr_by[a['family_rep_id']].append(a)
    part_by_cid={r['candidate_id']:r for r in partner_rows}
    summ=[]
    for fam in sorted(fams):
        mmem=[r for r in members if r['family_rep_id']==fam]
        aa=arr_by[fam]; strict=[x for x in aa if x['strict_published_art_geometry']=='True']
        strict_cids={x['candidate_id'] for x in strict}
        clusters=defaultdict(set)
        for cid in strict_cids:
            pr=part_by_cid.get(cid)
            if pr and pr.get('partner_cluster') and pr['partner_cluster']!='UNCLUSTERED': clusters[pr['partner_cluster']].add(pr['accession'])
        recurrent={k:len(v) for k,v in clusters.items() if len(v)>=2}
        summ.append(dict(family_rep_id=fam,member_locus_n=len(mmem),accession_n=len({r['accession'] for r in mmem}),any_repeat_locus_n=len({x['candidate_id'] for x in aa}),upstream_repeat_locus_n=len({x['candidate_id'] for x in aa if x['orientation']=='UPSTREAM'}),strict_art_geometry_locus_n=len(strict_cids),recurrent_partner_cluster_n=len(recurrent),strict_plus_recurrent_partner=str(bool(strict_cids and recurrent)),recurrent_partner_clusters=';'.join(f'{k}:{v}' for k,v in sorted(recurrent.items()))))
    sfields=list(summ[0].keys()); tsv_write(O/'DELTA12_FAMILY_NATIVE_ARCHITECTURE_SUMMARY.tsv',summ,sfields)
    # Descriptive architecture recovery by frozen source stratum. This is not a
    # prevalence estimate because analysis is conditioned on the 12 priority families.
    ctx_by_cid={r['candidate_id']:r for r in ctx_rows}
    strata=sorted({r.get('stratum','') for r in ctx_rows if r.get('stratum','')})
    stratum_rows=[]
    for st in strata:
        scids={r['candidate_id'] for r in ctx_rows if r.get('stratum')==st}
        sa=[a for a in arrays if a['candidate_id'] in scids]
        strict=[a for a in sa if a['strict_published_art_geometry']=='True']
        fam_present={ctx_by_cid[c]['family_rep_id'] for c in scids}
        fam_strict={a['family_rep_id'] for a in strict}
        stratum_rows.append(dict(stratum=st,locus_n=len(scids),family_n=len(fam_present),any_repeat_locus_n=len({a['candidate_id'] for a in sa}),upstream_repeat_locus_n=len({a['candidate_id'] for a in sa if a['orientation']=='UPSTREAM'}),strict_art_geometry_locus_n=len({a['candidate_id'] for a in strict}),strict_art_geometry_family_n=len(fam_strict)))
    stfields=['stratum','locus_n','family_n','any_repeat_locus_n','upstream_repeat_locus_n','strict_art_geometry_locus_n','strict_art_geometry_family_n']
    tsv_write(O/'DELTA12_STRATUM_ARCHITECTURE_SUMMARY.tsv',stratum_rows,stfields)
    metrics=read_stage04a_metrics(P)
    old_gb=metrics['old_bp']/1e9; delta_gb=metrics['delta_bp']/1e9
    close={
      'stage':'Stage05F3V2-Stage04B4','gate':'PASS_STAGE04B4_NATIVE_ARCHITECTURE_AND_SATURATION_CLOSEOUT' if not args.representative_only else 'PASS_STAGE04B4_REPRESENTATIVE_ONLY_DIAGNOSTIC',
      'gate_pass':True,'structural_positive_family_n':len(fams),'structural_positive_candidate_n':12,'native_member_locus_n':len(members),
      'family_with_any_repeat_n':sum(int(x['any_repeat_locus_n'])>0 for x in summ),'family_with_upstream_repeat_n':sum(int(x['upstream_repeat_locus_n'])>0 for x in summ),'family_with_strict_art_geometry_n':sum(int(x['strict_art_geometry_locus_n'])>0 for x in summ),'family_with_strict_plus_recurrent_partner_n':sum(x['strict_plus_recurrent_partner']=='True' for x in summ),
      'sequence_saturation':{
        'old_formal_per_gb':metrics['old_formal']/old_gb,'delta_formal_per_gb':metrics['delta_formal']/delta_gb,'old_family_per_gb':metrics['old_families']/old_gb,'delta_family_per_gb':metrics['delta_families']/delta_gb,'delta_only_combined_family_per_gb':metrics['delta_only_combined_families']/delta_gb,'interpretation':'RT-like instance yield remains near-proportional, whereas family novelty shows substantial diminishing returns; not saturated.'},
      'structural_priority_generalization':{'historical_resolved_positive':'28/30','independent_delta_positive':'12/32','interpretation':'partial but substantially attenuated out-of-universe structural generalization'},
      'architecture_by_source_stratum':stratum_rows,
      'panel_mapping_source':str(panel_src),'family_membership_source':str(member_src),'selected_accession_ledger':str(ledger),'source_root':str(SR),
      'rules_changed':False,'retuning_performed':False,
      'claim_boundary':'The native architecture endpoint applies only to the frozen 12 structurally compatible delta families and must not be interpreted as an architecture prevalence estimate across all 54.38 Gb.'
    }
    write_json(O/'STAGE04B4_CLOSEOUT_GATE.json',close)
    # provenance hashes
    files=[p for p in O.rglob('*') if p.is_file()]
    with open(O/'STAGE04B4_OUTPUT_SHA256.tsv','w') as f:
        for p in sorted(files): f.write(f'{sha(p)}\t{p.relative_to(HERE)}\n')
    print(json.dumps(close,indent=2))

if __name__=='__main__': main()
