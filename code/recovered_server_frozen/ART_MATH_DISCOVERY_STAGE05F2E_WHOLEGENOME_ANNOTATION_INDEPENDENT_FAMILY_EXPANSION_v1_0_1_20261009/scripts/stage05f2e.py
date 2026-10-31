#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math, os, re, shutil, subprocess, tarfile
from collections import defaultdict, Counter
from pathlib import Path

EXPECTED_QUERY_N=30
EXPECTED_GENOME_N=919
EXPECTED_TOTAL_BP=171397507
EXPECTED_SOURCE_MANIFEST_SHA="5eec6d7ee00584c6266b8a953b6112540148b949341e400099b285ce30d235f0"
EXPECTED_DUAL_SHA="04706102bdb8129e105c6cdba915d33e95e19963463869552622ad73283729ae"
MIN_AA=45
SEARCH_E=1e-2
CORE_E=1e-5
CORE_QCOV=0.70
CORE_FIDENT=0.25
BROAD_E=1e-3
BROAD_QCOV=0.50
MAX_SEQS=50000
FORMAL_E=1e-4
EXPLORE_E=1e-2
FLANK=10000

CODON_TABLE_11 = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
RC=str.maketrans('ACGTN','TGCAN')

F2A_NAMES=['STAGE05F2A_NOVELTY_FIRST_STRICT.faa','STAGE05F2A_NOVELTY_FIRST_STRICT_RESCUES.tsv','STAGE05F2A_NOVELTY_FIRST_INSTANCE_MAP.tsv','STAGE05F2A_GATE.json']
F2D_NAMES=['STAGE05F2D_ALL_SEARCH_HITS_ANNOTATED.tsv','STAGE05F2D_FAMILY_SUMMARY.tsv','STAGE05F2D_GATE.json']

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def sha_text(s): return hashlib.sha256(s.encode()).hexdigest()

def open_text(p):
    p=Path(p)
    return gzip.open(p,'rt') if p.suffix=='.gz' else p.open('r',newline='')

def read_tsv(p):
    with open_text(p) as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows,cols=None):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    if cols is None:
        cols=list(rows[0]) if rows else []
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore')
        if cols: w.writeheader()
        w.writerows(rows)

def asbool(x): return str(x).strip().lower() in ('1','true','yes','y')
def ff(x,default=float('nan')):
    try:return float(x)
    except:return default
def ii(x,default=0):
    try:return int(float(x))
    except:return default
def splitsemi(x): return [z for z in str(x or '').split(';') if z]

def resolve_dir_or_tar(x,names,work,tag):
    x=Path(x)
    if x.is_dir():
        out=x/'outputs'
        if all((out/n).exists() for n in names): return out,'live_directory'
    if x.is_file() and x.name.endswith(('.tar.gz','.tgz')):
        d=Path(work)/tag; d.mkdir(parents=True,exist_ok=True)
        wanted=set(names)
        with tarfile.open(x,'r:gz') as tf:
            for m in tf.getmembers():
                if not m.isfile(): continue
                if Path(m.name).name in wanted:
                    fh=tf.extractfile(m)
                    if fh: (d/Path(m.name).name).write_bytes(fh.read())
        if all((d/n).exists() for n in names): return d,'result_tar'
    raise RuntimeError(f'cannot resolve {names} from {x}')

def candidate_files(root,patterns):
    for pat in patterns:
        for p in Path(root).rglob(pat):
            if p.is_file(): yield p

def resolve_by_sha(root,expected,patterns,envkey):
    override=os.environ.get(envkey,'').strip()
    if override:
        p=Path(override)
        if sha256_file(p)!=expected: raise RuntimeError(f'{envkey} hash mismatch: {p}')
        return p,'override'
    seen=set()
    for p in candidate_files(root,patterns):
        rp=str(p.resolve())
        if rp in seen: continue
        seen.add(rp)
        try:
            if p.stat().st_size>250_000_000: continue
            if sha256_file(p)==expected: return p,'live_exact_sha256'
        except Exception: pass
    for p in sorted(Path(root).glob('**/*STAGE05B*.tar.gz'))+sorted(Path(root).glob('**/*STAGE05C1*.tar.gz')):
        try:
            with tarfile.open(p,'r:gz') as tf:
                for m in tf.getmembers():
                    if not m.isfile() or m.size>250_000_000: continue
                    bn=Path(m.name).name.lower()
                    if not any(x in bn for x in ('candidate','dual','manifest','region')): continue
                    fh=tf.extractfile(m)
                    if fh is None: continue
                    data=fh.read()
                    if hashlib.sha256(data).hexdigest()==expected:
                        d=Path(root)/'.stage05f2e_resolved'; d.mkdir(exist_ok=True)
                        out=d/Path(m.name).name; out.write_bytes(data)
                        return out,f'tar:{p.name}:{m.name}'
        except Exception: pass
    raise RuntimeError(f'could not resolve frozen file SHA256={expected}')

def parse_fasta(path,keep_ids=None):
    out={}; h=None; raw=''; seq=[]
    def emit():
        nonlocal h,raw,seq
        if h is None:return
        if keep_ids is None or h in keep_ids:
            s=''.join(seq).replace('*','').upper()
            out[h]={'seq':s,'sha256':sha_text(s),'aa_length':len(s),'header':raw}
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line:continue
            if line.startswith('>'):
                emit();raw=line[1:];h=raw.split()[0];seq=[]
            else:seq.append(line)
        emit()
    return out

def read_fasta_one(path):
    name=None;seq=[]
    with open_text(path) as f:
        for line in f:
            line=line.strip()
            if not line:continue
            if line.startswith('>'):
                if name is not None:break
                name=line[1:].split()[0]
            else:seq.append(line.upper())
    if name is None: raise RuntimeError(f'no FASTA record in {path}')
    return name,''.join(seq)

def resolve_fasta_manifest(root,manifest):
    rows=read_tsv(manifest); cols=list(rows[0]) if rows else []
    acc_col=next((x for x in ['sequence_id','accession','seq_id','genome_id'] if x in cols),None)
    path_col=next((x for x in cols if 'fasta' in x.lower() and ('path' in x.lower() or 'file' in x.lower() or 'frozen' in x.lower())),None)
    if not acc_col or not path_col: raise RuntimeError(f'cannot infer source manifest FASTA columns: {cols}')
    out={}
    for r in rows:
        acc=r[acc_col].strip(); q=r[path_col].strip()
        if not acc or not q:continue
        cands=[Path(q),Path(manifest).parent/q,Path(root)/q]
        cands+=list(Path(root).glob(f'**/{Path(q).name}'))[:4]
        hit=next((p for p in cands if p.exists() and p.is_file()),None)
        if hit:out[acc]=hit.resolve()
    return out,acc_col,path_col,rows

def detect_candidate_cols(rows):
    if not rows: raise RuntimeError('empty dual table')
    cols=list(rows[0])
    def choose(opts,req=True):
        for x in opts:
            if x in cols:return x
        if req:raise RuntimeError(f'missing {opts}; cols={cols}')
        return None
    return (choose(['sequence_id','accession','seq_id','genome_id']),
            choose(['start0','region_start0','union_start0','candidate_start0','start']),
            choose(['end0','region_end0','union_end0','candidate_end0','end']),
            choose(['region_id','candidate_id','union_region_id','id'],False))

def header(p):
    with open_text(p) as f:
        return next(csv.reader(f,delimiter='\t'),[])

def resolve_orf_registry(root):
    override=os.environ.get('STAGE05F2E_ORF_REGISTRY','').strip()
    if override and Path(override).exists():return Path(override),'override'
    cands=[]
    for b in Path(root).glob('ART_MATH_DISCOVERY_STAGE05C1*'):
        for p in b.rglob('*.tsv'):
            try:
                h=[x.lower() for x in header(p)]
                if any(x in h for x in ('sequence_id','accession')) and any(x in h for x in ('start0','cds_start0','orf_start0','protein_start0')) and any(x in h for x in ('end0','cds_end0','orf_end0','protein_end0')):
                    score=(10 if any(z in p.name.lower() for z in ('protein','orf','registry')) else 0)+min(p.stat().st_size/1e6,30)
                    cands.append((score,p))
            except:pass
    if not cands:return None,'not_resolved'
    cands.sort(reverse=True,key=lambda x:x[0])
    return cands[0][1],'auto'

def detect_orf_cols(rows):
    cols=list(rows[0]) if rows else []
    def choose(opts,req=True):
        for x in opts:
            if x in cols:return x
        if req:raise RuntimeError(f'missing ORF columns {opts}; have={cols}')
        return None
    return (choose(['sequence_id','accession','seq_id','genome_id']),
            choose(['start0','cds_start0','orf_start0','protein_start0','start']),
            choose(['end0','cds_end0','orf_end0','protein_end0','end']),
            choose(['strand','cds_strand','orf_strand'],False),
            choose(['record_key','protein_key','orf_id','target_name','id'],False))

def resolve_hmm(root):
    override=os.environ.get('STAGE05F2E_RVT_HMM','').strip()
    if override and Path(override).exists():return Path(override),'override'
    cands=[]
    for base in Path(root).glob('ART_MATH_DISCOVERY_STAGE05D1*'):
        for p in base.rglob('*.hmm'):
            try:
                n=p.read_text(errors='ignore').count('HMMER3/f')
                if 40<=n<=60:cands.append((abs(n-45),p,n))
            except:pass
    if not cands:
        for p in Path(root).rglob('*.hmm'):
            if 'RVT' not in p.name.upper() and 'RT' not in p.name.upper():continue
            try:
                n=p.read_text(errors='ignore').count('HMMER3/f')
                if 40<=n<=60:cands.append((abs(n-45),p,n))
            except:pass
    if not cands:raise RuntimeError('cannot resolve frozen ~45-profile RVT-All HMM; set STAGE05F2E_RVT_HMM')
    cands.sort(key=lambda x:(x[0],len(str(x[1]))))
    return cands[0][1],f'auto_models={cands[0][2]}'

def merge_intervals(items):
    out=[]
    for s,e in sorted(items):
        if not out or s>out[-1][1]:out.append([s,e])
        else:out[-1][1]=max(out[-1][1],e)
    return [(a,b) for a,b in out]

def overlaps_any(s,e,ints):
    # lists are short enough for hit-only post hoc use
    for a,b in ints:
        if b<=s:continue
        if a>=e:break
        if b>s and a<e:return True
    return False

def translate_frame(nt):
    return ''.join(CODON_TABLE_11.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))

def iter_whole_genome_segments(acc,seq):
    L=len(seq)
    for strand in ('+','-'):
        nt=seq if strand=='+' else seq.translate(RC)[::-1]
        for frame in (0,1,2):
            aa=translate_frame(nt[frame:])
            aa0=0
            for part in aa.split('*'):
                aa1=aa0+len(part)
                if len(part)>=MIN_AA:
                    n0=frame+aa0*3;n1=frame+aa1*3
                    if strand=='+':gs,ge=n0,n1
                    else:gs,ge=L-n1,L-n0
                    # post-hoc plausible start within first 30 aa
                    ntseg=seq[gs:ge] if strand=='+' else seq[gs:ge].translate(RC)[::-1]
                    first_start=-1
                    for j in range(0,min(len(part),30)):
                        cod=ntseg[j*3:j*3+3]
                        if cod in ('ATG','GTG','TTG'):
                            first_start=j;break
                    yield {'sequence_id':acc,'start0':gs,'end0':ge,'strand':strand,'frame':frame,
                           'aa_length':len(part),'aa_sequence':part,
                           'aa_sha256':sha_text(part),'genome_edge_partial':aa0==0 or aa1==len(aa),
                           'first_start_codon_aa_index':first_start,'start_codon_within_30aa':first_start>=0}
                aa0=aa1+1

def parse_search(path):
    fields=['query','target','evalue','bits','fident','alnlen','qcov','tcov','qlen','tlen']
    if not Path(path).exists() or Path(path).stat().st_size==0:return []
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t',fieldnames=fields))

def hit_class(r,exact):
    if exact:return 'EXACT_AA'
    e=ff(r['evalue'],9e99);q=ff(r['qcov'],0);fi=ff(r['fident'],0)
    if e<=CORE_E and q>=CORE_QCOV and fi>=CORE_FIDENT:return 'CORE_SEQUENCE_HOMOLOG'
    if e<=BROAD_E and q>=BROAD_QCOV:return 'BROAD_SEQUENCE_HOMOLOG_SENSITIVITY'
    return 'SEARCH_RESERVOIR'

def longer(ql,tl):return tl>=ql+30 and tl>=1.25*ql
def major(ql,tl):return longer(ql,tl) and tl>=max(180,1.5*ql)

def parse_domtbl(path,aa_id_to_sha):
    best={}
    if not Path(path).exists():return best
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):continue
            x=line.rstrip('\n').split(maxsplit=22)
            if len(x)<23:continue
            target=x[0];tlen=ii(x[2]);query=x[3];qlen=ii(x[5])
            full_e=ff(x[6],9e99);dom_i=ff(x[12],9e99);score=ff(x[13],float('-inf'))
            hmm_from=ii(x[15]);hmm_to=ii(x[16]);ali_from=ii(x[17]);ali_to=ii(x[18]);acc=ff(x[21],float('nan'))
            sha=aa_id_to_sha.get(query)
            if not sha:continue
            row={'best_rvt_target':target,'best_full_evalue':full_e,'best_domain_i_evalue':dom_i,'domain_score':score,
                 'query_coverage':((ali_to-ali_from+1)/qlen if qlen else None),
                 'model_coverage':((hmm_to-hmm_from+1)/tlen if tlen else None),
                 'domain_acc':acc,'hmm_from':hmm_from,'hmm_to':hmm_to,'ali_from':ali_from,'ali_to':ali_to}
            key=(full_e,dom_i,-score)
            if sha not in best or key<best[sha][0]:best[sha]=(key,row)
    return {k:v[1] for k,v in best.items()}

class UF:
    def __init__(self):self.p={}
    def add(self,x):
        if x not in self.p:self.p[x]=x
    def find(self,x):
        self.add(x)
        while self.p[x]!=x:
            self.p[x]=self.p[self.p[x]];x=self.p[x]
        return x
    def union(self,a,b):
        ra,rb=self.find(a),self.find(b)
        if ra!=rb:self.p[rb]=ra

def build_f2d_components(f2d_hits,query_ids):
    uf=UF()
    for q in query_ids:uf.add(q)
    bysha=defaultdict(set)
    for r in f2d_hits:
        if r.get('hit_class')!='CORE_SEQUENCE_HOMOLOG':continue
        q=r.get('query','');t=r.get('target','');h=r.get('target_aa_sha256','')
        if q in query_ids:
            bysha[h].add(q)
            if t in query_ids:uf.union(q,t)
    for qs in bysha.values():
        qs=sorted(qs)
        for q in qs[1:]:uf.union(qs[0],q)
    groups=defaultdict(list)
    for q in sorted(query_ids):groups[uf.find(q)].append(q)
    groups_sorted=sorted(groups.values(),key=lambda x:(x[0],len(x)))
    out={}
    for i,g in enumerate(groups_sorted,1):
        cid=f'F2D_COMP_{i:03d}'
        for q in g:out[q]=cid
    return out,groups_sorted

def run(args):
    root=Path(args.project_root).resolve();work=Path(args.work);out=Path(args.outdir);raw=Path(args.rawdir);logs=Path(args.logs)
    for d in (work,out,raw,logs):d.mkdir(parents=True,exist_ok=True)
    if not shutil.which('mmseqs'):raise RuntimeError('mmseqs not found in PATH')
    if not shutil.which('hmmscan'):raise RuntimeError('hmmscan not found in PATH')
    f2a,f2a_route=resolve_dir_or_tar(args.f2a,F2A_NAMES,work,'f2a')
    f2d,f2d_route=resolve_dir_or_tar(args.f2d,F2D_NAMES,work,'f2d')
    g2a=json.loads((f2a/'STAGE05F2A_GATE.json').read_text());g2d=json.loads((f2d/'STAGE05F2D_GATE.json').read_text())
    if not g2a.get('gate_pass') or not g2d.get('gate_pass'):raise RuntimeError('upstream F2A/F2D gate not PASS')
    qrows=read_tsv(f2a/'STAGE05F2A_NOVELTY_FIRST_STRICT_RESCUES.tsv');qseq=parse_fasta(f2a/'STAGE05F2A_NOVELTY_FIRST_STRICT.faa')
    qinst=read_tsv(f2a/'STAGE05F2A_NOVELTY_FIRST_INSTANCE_MAP.tsv')
    if len(qrows)!=EXPECTED_QUERY_N or len(qseq)!=EXPECTED_QUERY_N:raise RuntimeError('F2A query count !=30')
    qmeta={r['query_id']:r for r in qrows}
    if set(qmeta)!=set(qseq):raise RuntimeError('F2A query IDs mismatch')
    for q in qseq:
        if qseq[q]['sha256']!=qmeta[q].get('aa_sha256',''):raise RuntimeError(f'F2A hash mismatch {q}')
    qhash={q:x['sha256'] for q,x in qseq.items()}
    qids_by_hash=defaultdict(list)
    for q,h in qhash.items():qids_by_hash[h].append(q)

    mani,mani_route=resolve_by_sha(root,EXPECTED_SOURCE_MANIFEST_SHA,['*MANIFEST*.tsv','*manifest*.tsv'],'STAGE05F2E_SOURCE_MANIFEST')
    dual,dual_route=resolve_by_sha(root,EXPECTED_DUAL_SHA,['*DUAL*.tsv','*dual*.tsv','*CANDIDATE*.tsv','*candidate*.tsv','*REGION*.tsv'],'STAGE05F2E_DUAL_TSV')
    fmap,mac,mpath,mrows=resolve_fasta_manifest(root,mani)
    if len(mrows)!=EXPECTED_GENOME_N:raise RuntimeError(f'source manifest rows={len(mrows)}, expected919')
    accs=[r[mac].strip() for r in mrows if r[mac].strip()]
    miss=[a for a in accs if a not in fmap]
    if miss:raise RuntimeError(f'missing {len(miss)} frozen FASTAs; first={miss[:10]}')
    hmm,hmm_route=resolve_hmm(root)
    orf,orf_route=resolve_orf_registry(root)

    # v1.0.1 hotfix preflight: replay the frozen F2A instance coordinates directly
    # against the same frozen 919-genome source FASTAs. This is independent of MMseqs.
    # A strict non-edge F2A peptide must be exactly reconstructable from at least one
    # listed genomic instance before any whole-genome family search is interpreted.
    coord_ok=set();coord_bad=[];coord_instance_n=0
    src_cache={}
    for r in qinst:
        q=r.get('query_id','')
        if q not in qseq:continue
        a=r.get('sequence_id','')
        if a not in fmap:
            coord_bad.append({'query_id':q,'sequence_id':a,'reason':'ACCESSION_NOT_IN_FROZEN_MANIFEST'});continue
        if a not in src_cache:src_cache[a]=read_fasta_one(fmap[a])[1]
        g=src_cache[a];st=ii(r.get('start0'));en=ii(r.get('end0'));strand=r.get('strand','')
        nt=g[st:en]
        if strand=='-':nt=nt.translate(RC)[::-1]
        aa=translate_frame(nt)
        coord_instance_n+=1
        if aa==qseq[q]['seq']:
            coord_ok.add(q)
        else:
            coord_bad.append({'query_id':q,'sequence_id':a,'start0':st,'end0':en,'strand':strand,
                              'expected_sha256':qhash[q],'replayed_sha256':sha_text(aa),'reason':'COORDINATE_REPLAY_AA_MISMATCH'})
    if set(coord_ok)!=set(qseq):
        missing=sorted(set(qseq)-coord_ok)
        raise RuntimeError(f'F2A coordinate replay failed for {len(missing)}/{EXPECTED_QUERY_N} queries; first={missing[:10]}')
    if coord_bad:
        raise RuntimeError(f'F2A instance coordinate replay had {len(coord_bad)} mismatching instances; first={coord_bad[:3]}')

    # Frozen candidate +/-10kb windows for outside-neighborhood flag.
    drows=read_tsv(dual);dac,dst,den,drid=detect_candidate_cols(drows)
    if len(drows)!=3674:raise RuntimeError(f'dual rows={len(drows)}, expected3674')
    win_by_acc=defaultdict(list)
    genome_len_cache={}
    for a in sorted(set(r[dac] for r in drows)):
        _,s=read_fasta_one(fmap[a]);genome_len_cache[a]=len(s)
    for r in drows:
        a=r[dac];s=ii(r[dst]);e=ii(r[den]);L=genome_len_cache[a]
        win_by_acc[a].append((max(0,s-FLANK),min(L,e+FLANK)))
    for a in list(win_by_acc):win_by_acc[a]=merge_intervals(win_by_acc[a])

    # Build ORF map post hoc.
    orf_by_acc=defaultdict(list)
    if orf:
        orows=read_tsv(orf)
        oa,os_,oe,ostr,okey=detect_orf_cols(orows)
        for r in orows:
            try:orf_by_acc[r[oa]].append((ii(r[os_]),ii(r[oe]),r.get(ostr,'') if ostr else '',r.get(okey,'') if okey else ''))
            except:pass
        for a in orf_by_acc:orf_by_acc[a].sort()

    # F2D query components + semantic priority audit.
    f2dh=read_tsv(f2d/'STAGE05F2D_ALL_SEARCH_HITS_ANNOTATED.tsv')
    comp_map,groups=build_f2d_components(f2dh,set(qseq))
    f2df=read_tsv(f2d/'STAGE05F2D_FAMILY_SUMMARY.tsv')
    label_audit=[]
    for r in f2df:
        core=ii(r.get('core_homolog_unique_aa_n'));coreg=ii(r.get('core_homolog_genome_n'));maj=ii(r.get('major_extension_unique_aa_n'));formalmaj=ii(r.get('formal_rvt_major_extension_unique_aa_n'))
        if formalmaj and ii(r.get('major_extension_genome_n'))>=2:new='P1_LONGER_FORMAL_RVT_RESCUE_MULTI_GENOME'
        elif formalmaj:new='P2_LONGER_FORMAL_RVT_RESCUE_SINGLE_GENOME'
        elif maj:new='P3_MAJOR_EXTENSION_NO_FORMAL_RVT'
        elif core and coreg>=2:new='P4A_REPLICATED_CORE_HOMOLOG_NO_MAJOR_EXTENSION'
        elif core:new='P4B_SINGLE_CORE_HOMOLOG_NO_MAJOR_EXTENSION'
        else:new='P5_NO_CORE_HOMOLOG'
        label_audit.append({'query_id':r['query_id'],'original_f2d_priority':r.get('f2d_priority',''),'semantic_r1_priority':new,
                            'core_homolog_unique_aa_n':core,'core_homolog_genome_n':coreg,'major_extension_unique_aa_n':maj})
    write_tsv(out/'STAGE05F2E_F2D_PRIORITY_LABEL_AUDIT.tsv',label_audit)

    # Whole-genome annotation-independent peptide generation.
    # v1.0.1 can reuse the completed v1.0.0 whole-genome peptide materialization,
    # because v1.0.0 failed only after generation/search at the exact-recovery check.
    wholefa=work/'STAGE05F2E_WHOLE919_SIXFRAME_INSTANCES.faa'
    inst=work/'STAGE05F2E_WHOLE919_SIXFRAME_INSTANCES.tsv'
    total_bp=0
    for a in accs:
        _,seq=read_fasta_one(fmap[a]);total_bp+=len(seq)
    if len(accs)!=EXPECTED_GENOME_N or total_bp!=EXPECTED_TOTAL_BP:raise RuntimeError(f'genome universe mismatch n={len(accs)} bp={total_bp}')

    reuse=False
    if args.reuse_v100_work:
        rw=Path(args.reuse_v100_work)
        rfa=rw/'STAGE05F2E_WHOLE919_SIXFRAME_INSTANCES.faa'
        rit=rw/'STAGE05F2E_WHOLE919_SIXFRAME_INSTANCES.tsv'
        if rfa.exists() and rfa.stat().st_size>0 and rit.exists() and rit.stat().st_size>0:
            for src,dst in ((rfa,wholefa),(rit,inst)):
                if dst.exists():dst.unlink()
                try:os.link(src,dst)
                except OSError:shutil.copy2(src,dst)
            reuse=True
    if not reuse:
        pepn=0
        with wholefa.open('w') as ffh, inst.open('w',newline='') as tfh:
            cols=['peptide_instance_id','sequence_id','start0','end0','strand','frame','aa_length','aa_sha256','genome_edge_partial','first_start_codon_aa_index','start_codon_within_30aa']
            w=csv.DictWriter(tfh,fieldnames=cols,delimiter='\t');w.writeheader()
            for a in accs:
                _,seq=read_fasta_one(fmap[a])
                for p in iter_whole_genome_segments(a,seq):
                    pepn+=1;pid=f'WG6F_{pepn:09d}'
                    ffh.write(f'>{pid} acc={a} start0={p["start0"]} end0={p["end0"]} strand={p["strand"]} frame={p["frame"]} sha256={p["aa_sha256"]}\n{p["aa_sequence"]}\n')
                    rr={k:p.get(k,'') for k in cols};rr['peptide_instance_id']=pid;w.writerow(rr)

    # Direct exact-AA recovery audit from the materialized whole-genome instance registry.
    # This is the correct invariant. It does NOT depend on a heuristic similarity-search prefilter.
    direct_exact=set();exact_target_ids=defaultdict(list);pepn=0
    with inst.open(newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            pepn+=1;h=r.get('aa_sha256','')
            for q in qids_by_hash.get(h,[]):
                direct_exact.add(q);exact_target_ids[q].append(r['peptide_instance_id'])
    if set(direct_exact)!=set(qseq):
        missing=sorted(set(qseq)-direct_exact)
        raise RuntimeError(f'direct whole-genome AA-hash recovery {len(direct_exact)}/{EXPECTED_QUERY_N}; missing={missing[:15]}')

    # MMseqs family search. v1.0.0 incorrectly used MMseqs-returned exact hits as the
    # generation invariant. For short/low-complexity phage peptides, default masking can
    # suppress even an exact target from the returned hit list. v1.0.1 disables masking
    # for this frozen 30-query search and uses an exhaustive fallback only for queries
    # whose directly verified exact target is still absent from the primary result.
    qfa=f2a/'STAGE05F2A_NOVELTY_FIRST_STRICT.faa'
    primary=raw/'STAGE05F2E_QUERY30_VS_WHOLE919_SIXFRAME_PRIMARY_MASK0.m8'
    m8=raw/'STAGE05F2E_QUERY30_VS_WHOLE919_SIXFRAME.m8'
    tmp=work/'mmseqs_tmp_primary'
    if tmp.exists():shutil.rmtree(tmp)
    cmd=['mmseqs','easy-search',str(qfa),str(wholefa),str(primary),str(tmp),'-s','7.5','-e',str(SEARCH_E),'--max-seqs',str(MAX_SEQS),'--mask','0','--threads',str(args.threads),
         '--format-output','query,target,evalue,bits,fident,alnlen,qcov,tcov,qlen,tlen']
    with (logs/'mmseqs_wholegenome_primary_mask0.log').open('w') as lf:
        subprocess.run(cmd,check=True,stdout=lf,stderr=subprocess.STDOUT)
    hs_primary=parse_search(primary)
    primary_targets=defaultdict(set)
    for r in hs_primary:primary_targets[r['query']].add(r['target'])
    mm_exact_primary={q for q in qseq if any(t in primary_targets.get(q,set()) for t in exact_target_ids[q])}
    missing_primary=sorted(set(qseq)-mm_exact_primary)

    hs_fallback=[]
    if missing_primary:
        mfa=work/'STAGE05F2E_MISSING_EXACT_PRIMARY_QUERIES.faa'
        with mfa.open('w') as f:
            for q in missing_primary:f.write(f'>{q}\n{qseq[q]["seq"]}\n')
        fb=raw/'STAGE05F2E_QUERY_MISSING_PRIMARY_EXHAUSTIVE.m8'
        tmp2=work/'mmseqs_tmp_exhaustive'
        if tmp2.exists():shutil.rmtree(tmp2)
        cmd=['mmseqs','easy-search',str(mfa),str(wholefa),str(fb),str(tmp2),'-e',str(SEARCH_E),'--max-seqs',str(MAX_SEQS),'--mask','0','--exhaustive-search','1','--threads',str(args.threads),
             '--format-output','query,target,evalue,bits,fident,alnlen,qcov,tcov,qlen,tlen']
        with (logs/'mmseqs_wholegenome_exhaustive_fallback.log').open('w') as lf:
            subprocess.run(cmd,check=True,stdout=lf,stderr=subprocess.STDOUT)
        hs_fallback=parse_search(fb)

    # Deterministically merge primary and fallback rows.
    fields=['query','target','evalue','bits','fident','alnlen','qcov','tcov','qlen','tlen']
    merged=[];seen=set()
    for r in hs_primary+hs_fallback:
        k=tuple(r.get(x,'') for x in fields)
        if k in seen:continue
        seen.add(k);merged.append(r)
    with m8.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writerows(merged)
    hs=merged
    byq=defaultdict(list);target_ids=set()
    for r in hs:byq[r['query']].append(r);target_ids.add(r['target'])
    search_cap=any(len(byq[q])>=MAX_SEQS for q in qseq)
    mm_exact_final={q for q in qseq if any(t in set(x['target'] for x in byq.get(q,[])) for t in exact_target_ids[q])}
    if set(mm_exact_final)!=set(qseq):
        missing=sorted(set(qseq)-mm_exact_final)
        raise RuntimeError(f'MMseqs exact-target recovery remained incomplete after mask-off/exhaustive fallback: {len(mm_exact_final)}/{EXPECTED_QUERY_N}; missing={missing[:15]}')

    # Recover metadata and sequences only for returned target instances.
    meta={}
    with inst.open(newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['peptide_instance_id'] in target_ids:meta[r['peptide_instance_id']]=r
    if set(meta)!=target_ids:raise RuntimeError(f'missing hit metadata {len(target_ids-set(meta))}')
    seqs={}
    h=None;seq=[]
    with wholefa.open() as f:
        for line in f:
            line=line.strip()
            if not line:continue
            if line.startswith('>'):
                if h in target_ids:seqs[h]=''.join(seq)
                h=line[1:].split()[0];seq=[]
            else:
                if h in target_ids:seq.append(line)
        if h in target_ids:seqs[h]=''.join(seq)
    if set(seqs)!=target_ids:raise RuntimeError(f'missing hit sequences {len(target_ids-set(seqs))}')

    # Final exact-AA audit: direct generation and MMseqs are recorded separately.
    exact_rec=set()
    for q,rs in byq.items():
        for r in rs:
            if sha_text(seqs[r['target']])==qhash[q]:exact_rec.add(q)
    if set(exact_rec)!=set(qseq):raise RuntimeError(f'final exact MMseqs sequence recovery {len(exact_rec)}/{EXPECTED_QUERY_N}')
    audit=[]
    for q in sorted(qseq):
        audit.append({'query_id':q,'aa_sha256':qhash[q],'coordinate_replay_pass':str(q in coord_ok),
                      'direct_wholegenome_hash_recovery':str(q in direct_exact),'direct_exact_instance_n':len(exact_target_ids[q]),
                      'mmseqs_primary_mask0_exact_recovery':str(q in mm_exact_primary),
                      'mmseqs_final_exact_recovery':str(q in mm_exact_final),
                      'exhaustive_fallback_used':str(q in missing_primary)})
    write_tsv(out/'STAGE05F2E_EXACT_QUERY_RECOVERY_AUDIT.tsv',audit)

    # Deduplicate returned target AAs and HMM-scan them.
    sha_to_seq={};sha_to_aid={};aid_to_sha={}
    for t,s in seqs.items():
        hsh=sha_text(s)
        if hsh not in sha_to_seq:
            aid=f'WGAA_{len(sha_to_seq)+1:07d}'
            sha_to_seq[hsh]=s;sha_to_aid[hsh]=aid;aid_to_sha[aid]=hsh
    hitfa=out/'STAGE05F2E_RETURNED_TARGET_UNIQUE_AA.faa'
    with hitfa.open('w') as f:
        for hsh in sorted(sha_to_seq):
            aid=sha_to_aid[hsh];f.write(f'>{aid} sha256={hsh} len={len(sha_to_seq[hsh])}\n{sha_to_seq[hsh]}\n')
    dom=raw/'STAGE05F2E_RETURNED_TARGET_RVT_ALL.domtblout';tbl=raw/'STAGE05F2E_RETURNED_TARGET_RVT_ALL.tblout'
    cmd=['hmmscan','--cpu',str(args.threads),'-E',str(EXPLORE_E),'--domE',str(EXPLORE_E),'--domtblout',str(dom),'--tblout',str(tbl),str(hmm),str(hitfa)]
    with (logs/'hmmscan_returned_targets.log').open('w') as lf:
        subprocess.run(cmd,check=True,stdout=lf,stderr=subprocess.STDOUT)
    hmm_best=parse_domtbl(dom,aid_to_sha)

    # Helper: Prodigal relationship only for hit instances.
    def prod_rel(r):
        if not orf:return 'ORF_REGISTRY_UNRESOLVED',0
        a=r['sequence_id'];s=ii(r['start0']);e=ii(r['end0']);strand=r['strand'];same=[]
        for x,y,st,key in orf_by_acc.get(a,[]):
            if y<=s:continue
            if x>=e:break
            if y>s and x<e and (not st or st==strand):same.append((x,y,st,key))
        if not same:return 'NO_SAME_STRAND_PRODIGAL_OVERLAP',0
        if any(x<=s and y>=e for x,y,_,_ in same):return 'CONTAINED_IN_EXISTING_PRODIGAL_ORF',len(same)
        return 'BOUNDARY_OR_PARTIAL_DISCORDANT_WITH_PRODIGAL',len(same)

    annot=[]
    for q in qseq:
        for r in byq.get(q,[]):
            t=r['target'];m=meta[t];s=seqs[t];hsh=sha_text(s);exact=(hsh==qhash[q]);hc=hit_class(r,exact)
            ql=ii(r['qlen']);tl=ii(r['tlen']);lg=(not exact and hc=='CORE_SEQUENCE_HOMOLOG' and longer(ql,tl));mj=lg and major(ql,tl)
            hb=hmm_best.get(hsh,{})
            formal=bool(hb and hb.get('best_full_evalue',9e99)<=FORMAL_E and hb.get('best_domain_i_evalue',9e99)<=FORMAL_E)
            explore=bool(hb and hb.get('best_full_evalue',9e99)<=EXPLORE_E and hb.get('best_domain_i_evalue',9e99)<=EXPLORE_E)
            pr,nolap=prod_rel(m)
            inside=overlaps_any(ii(m['start0']),ii(m['end0']),win_by_acc.get(m['sequence_id'],[]))
            x=dict(r);x.update({
                'query_aa_sha256':qhash[q],'target_aa_sha256':hsh,'hit_class':hc,'exact_aa':str(exact),
                'longer_core_homolog':str(lg),'major_extension':str(mj),
                'target_sequence_id':m['sequence_id'],'target_start0':m['start0'],'target_end0':m['end0'],'target_strand':m['strand'],'target_frame':m['frame'],
                'target_aa_length':m['aa_length'],'target_genome_edge_partial':m['genome_edge_partial'],
                'target_start_codon_within_30aa':m['start_codon_within_30aa'],'target_first_start_codon_aa_index':m['first_start_codon_aa_index'],
                'target_inside_frozen_f1_neighborhood':str(inside),'target_outside_frozen_f1_neighborhood':str(not inside),
                'target_prodigal_relation':pr,'target_same_strand_prodigal_overlap_n':nolap,
                'target_formal_rvt_profile':str(formal),'target_exploratory_rvt_profile':str(explore),
                'target_best_rvt_target':hb.get('best_rvt_target',''),'target_best_full_evalue':hb.get('best_full_evalue',''),
                'target_best_domain_i_evalue':hb.get('best_domain_i_evalue',''),'target_query_coverage':hb.get('query_coverage',''),
                'target_model_coverage':hb.get('model_coverage',''),'target_domain_score':hb.get('domain_score','')
            })
            annot.append(x)
    write_tsv(out/'STAGE05F2E_ALL_WHOLEGENOME_HITS_ANNOTATED.tsv',annot)

    # Per-query summaries, deduplicating by AA for sequence counts and accession for genome counts.
    summaries=[];longer_rows=[];formal_longer_rows=[]
    for q in sorted(qseq):
        rs=[x for x in annot if x['query']==q]
        core=[x for x in rs if x['hit_class']=='CORE_SEQUENCE_HOMOLOG']
        broad=[x for x in rs if x['hit_class']=='BROAD_SEQUENCE_HOMOLOG_SENSITIVITY']
        ex=[x for x in core if asbool(x['longer_core_homolog'])]
        fex=[x for x in ex if asbool(x['target_formal_rvt_profile'])]
        out_fex=[x for x in fex if asbool(x['target_outside_frozen_f1_neighborhood'])]
        maj=[x for x in ex if asbool(x['major_extension'])]
        fmaj=[x for x in maj if asbool(x['target_formal_rvt_profile'])]
        def genomes(rows):return set(x['target_sequence_id'] for x in rows)
        def uniqaa(rows):return set(x['target_aa_sha256'] for x in rows)
        # preserve every instance for downstream audit
        for x in ex:
            y=dict(x);y['f2d_component_id']=comp_map.get(q,'');longer_rows.append(y)
            if asbool(x['target_formal_rvt_profile']):formal_longer_rows.append(y)
        longest=max(core,key=lambda x:ii(x['target_aa_length'])) if core else None
        if out_fex:
            pri='Q1_OUTSIDE_F1_NEIGHBORHOOD_LONGER_FORMAL_RVT'
        elif fex:
            pri='Q2_WITHIN_F1_NEIGHBORHOOD_LONGER_FORMAL_RVT'
        elif any(asbool(x['target_outside_frozen_f1_neighborhood']) for x in ex):
            pri='Q3_OUTSIDE_F1_NEIGHBORHOOD_LONGER_HMM_NEGATIVE'
        elif ex:
            pri='Q4_WITHIN_NEIGHBORHOOD_LONGER_HMM_NEGATIVE'
        elif core:
            pri='Q5_CORE_FAMILY_NO_LONGER_HOMOLOG'
        else:
            pri='Q6_NO_CORE_WHOLEGENOME_EXPANSION'
        summaries.append({
            'query_id':q,'f2d_component_id':comp_map.get(q,''),'query_aa_length':qseq[q]['aa_length'],
            'query_best_rvt_target':qmeta[q].get('best_rvt_target',''),'query_best_full_evalue':qmeta[q].get('best_full_evalue',''),
            'query_query_coverage':qmeta[q].get('query_coverage',''),'query_model_coverage':qmeta[q].get('model_coverage',''),
            'query_yxdd_like':qmeta[q].get('yxdd_like',''),'wholegenome_search_hit_n':len(rs),
            'core_homolog_unique_aa_n':len(uniqaa(core)),'core_homolog_genome_n':len(genomes(core)),
            'broad_homolog_unique_aa_n':len(uniqaa(broad)),'broad_homolog_genome_n':len(genomes(broad)),
            'longer_core_unique_aa_n':len(uniqaa(ex)),'longer_core_genome_n':len(genomes(ex)),
            'formal_longer_unique_aa_n':len(uniqaa(fex)),'formal_longer_genome_n':len(genomes(fex)),
            'outside_f1_formal_longer_unique_aa_n':len(uniqaa(out_fex)),'outside_f1_formal_longer_genome_n':len(genomes(out_fex)),
            'major_extension_unique_aa_n':len(uniqaa(maj)),'formal_major_extension_unique_aa_n':len(uniqaa(fmaj)),
            'longest_core_target':longest['target'] if longest else '','longest_core_target_len':longest['target_aa_length'] if longest else '',
            'f2e_priority':pri
        })
    write_tsv(out/'STAGE05F2E_FAMILY_SUMMARY.tsv',summaries)
    write_tsv(out/'STAGE05F2E_LONGER_WHOLEGENOME_HOMOLOG_INSTANCES.tsv',longer_rows)
    write_tsv(out/'STAGE05F2E_FORMAL_LONGER_WHOLEGENOME_HOMOLOG_INSTANCES.tsv',formal_longer_rows)

    # Component summaries.
    comp_rows=[]
    comp_ids=sorted(set(comp_map.values()))
    for cid in comp_ids:
        qs=sorted(q for q,c in comp_map.items() if c==cid)
        ss=[x for x in summaries if x['query_id'] in qs]
        ann=[x for x in annot if x['query'] in qs]
        ex=[x for x in ann if asbool(x['longer_core_homolog'])]
        fex=[x for x in ex if asbool(x['target_formal_rvt_profile'])]
        of=[x for x in fex if asbool(x['target_outside_frozen_f1_neighborhood'])]
        maj=[x for x in ex if asbool(x['major_extension'])]
        fmaj=[x for x in maj if asbool(x['target_formal_rvt_profile'])]
        comp_rows.append({
            'f2d_component_id':cid,'query_n':len(qs),'query_ids':';'.join(qs),
            'query_profiles':';'.join(sorted(set(qmeta[q].get('best_rvt_target','') for q in qs if qmeta[q].get('best_rvt_target','')))),
            'wholegenome_core_unique_aa_n':len(set(x['target_aa_sha256'] for x in ann if x['hit_class']=='CORE_SEQUENCE_HOMOLOG')),
            'wholegenome_core_genome_n':len(set(x['target_sequence_id'] for x in ann if x['hit_class']=='CORE_SEQUENCE_HOMOLOG')),
            'longer_unique_aa_n':len(set(x['target_aa_sha256'] for x in ex)),
            'formal_longer_unique_aa_n':len(set(x['target_aa_sha256'] for x in fex)),
            'outside_f1_formal_longer_unique_aa_n':len(set(x['target_aa_sha256'] for x in of)),
            'major_extension_unique_aa_n':len(set(x['target_aa_sha256'] for x in maj)),
            'formal_major_extension_unique_aa_n':len(set(x['target_aa_sha256'] for x in fmaj))
        })
    write_tsv(out/'STAGE05F2E_QUERY_COMPONENT_SUMMARY.tsv',comp_rows)

    # FASTA of unique formal longer targets.
    chosen={}
    for x in formal_longer_rows:
        chosen[x['target_aa_sha256']]=seqs[x['target']]
    with (out/'STAGE05F2E_FORMAL_LONGER_TARGETS.faa').open('w') as f:
        for hsh in sorted(chosen):
            f.write(f'>F2E_FORMAL_LONGER_{hsh[:16]} sha256={hsh} len={len(chosen[hsh])}\n{chosen[hsh]}\n')

    semcnt=Counter(x['semantic_r1_priority'] for x in label_audit)
    prcnt=Counter(x['f2e_priority'] for x in summaries)
    component_primary=sum(ii(x['outside_f1_formal_longer_unique_aa_n'])>0 for x in comp_rows)
    gate={
        'stage':'Stage05F2E','version':'v1.0.1','gate':'PASS_STAGE05F2E_WHOLEGENOME_ANNOTATION_INDEPENDENT_FAMILY_EXPANSION','gate_pass':True,
        'f2a_route':f2a_route,'f2d_route':f2d_route,'source_manifest_route':mani_route,'dual_candidate_route':dual_route,
        'source_manifest_sha256':sha256_file(mani),'dual_candidate_sha256':sha256_file(dual),
        'genome_n':len(accs),'total_bp':total_bp,'wholegenome_sixframe_instance_n':pepn,
        'query_n':len(qseq),'f2d_query_component_n':len(comp_ids),'mmseqs_hit_rows_n':len(hs),'mmseqs_search_cap_reached':search_cap,
        'coordinate_replay_query_n':len(coord_ok),'coordinate_replay_instance_n':coord_instance_n,'direct_wholegenome_hash_recovery_n':len(direct_exact),
        'mmseqs_primary_mask0_exact_recovery_n':len(mm_exact_primary),'mmseqs_exhaustive_fallback_query_n':len(missing_primary),
        'mmseqs_final_exact_recovery_n':len(mm_exact_final),'exact_query_recovery_n':len(exact_rec),'wholegenome_materialization_reused_from_v100':reuse,'returned_target_unique_aa_n':len(sha_to_seq),
        'wholegenome_core_expansion_query_n':sum(ii(x['core_homolog_unique_aa_n'])>0 for x in summaries),
        'longer_homolog_query_n':sum(ii(x['longer_core_unique_aa_n'])>0 for x in summaries),
        'formal_longer_homolog_query_n':sum(ii(x['formal_longer_unique_aa_n'])>0 for x in summaries),
        'outside_f1_neighborhood_formal_longer_query_n':sum(ii(x['outside_f1_formal_longer_unique_aa_n'])>0 for x in summaries),
        'outside_f1_neighborhood_formal_longer_component_n':component_primary,
        'formal_major_extension_query_n':sum(ii(x['formal_major_extension_unique_aa_n'])>0 for x in summaries),
        'f2d_semantic_label_audit_counts':dict(semcnt),'f2e_priority_counts':dict(prcnt),
        'new_hmm_search_performed':True,'new_structure_prediction_performed':False,'implementation_hotfix':'v1.0.0 exact-recovery invariant moved from heuristic MMseqs return to direct whole-genome AA-hash audit; MMseqs masking disabled with exhaustive fallback for exact-target completeness',
        'annotation_used_to_generate_peptides':False,'annotation_used_post_hoc_only':bool(orf),
        'orf_registry_route':orf_route,'rvt_hmm_route':hmm_route,
        'dna_score_recomputed':False,'dna_rank_recomputed':False,'candidate_membership_changed':False,
        'claim_boundary':'Whole-genome annotation-independent family expansion plus frozen-threshold HMM triage of returned homologs. No result establishes biochemical RT activity or ART identity; no-hit does not exclude more divergent/recoded systems.'
    }
    (out/'STAGE05F2E_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    notes=[
      '# Stage05F2E manuscript-facing notes\n\n',
      f'- Frozen source universe: {len(accs)} genomes / {total_bp:,} bp.\n',
      f'- Whole-genome six-frame stop-to-stop peptide instances (>=45 aa): {pepn:,}.\n',
      f'- Frozen novelty-first queries: {len(qseq)}; coordinate replay: {len(coord_ok)}/{len(qseq)}; direct whole-genome AA-hash recovery: {len(direct_exact)}/{len(qseq)}; final MMseqs exact-target recovery: {len(mm_exact_final)}/{len(qseq)}.\n',
      f'- F2D query components after collapsing connected sequence families: {len(comp_ids)}.\n',
      f'- Queries with longer core homologs anywhere in the 919 genomes: {gate["longer_homolog_query_n"]}.\n',
      f'- Queries with formal-RVT longer homologs: {gate["formal_longer_homolog_query_n"]}.\n',
      f'- Query-components with formal-RVT longer homologs outside the frozen F1 candidate-neighborhood universe: {component_primary}.\n',
      f'- Queries with formal-RVT major extensions: {gate["formal_major_extension_query_n"]}.\n\n',
      'Primary interpretation must use component-level counts because several F2A queries are members of the same sequence family. '
      'Prodigal relations are post-hoc only. Whole-genome longer homologs remain representation evidence, not biochemical RT calls.\n'
    ]
    (out/'STAGE05F2E_MANUSCRIPT_NOTES.md').write_text(''.join(notes))

    # checksums of compact outputs; omit huge work FASTA/TSV.
    files=[m8,dom,tbl,hitfa]
    files+=sorted(out.glob('STAGE05F2E_*.tsv'))+sorted(out.glob('STAGE05F2E_*.json'))+sorted(out.glob('STAGE05F2E_*.md'))+sorted(out.glob('STAGE05F2E_*.faa'))
    seen=set();rows=[]
    for p in files:
        p=Path(p)
        if p.name=='STAGE05F2E_SHA256.tsv':continue
        if not p.exists() or p in seen:continue
        seen.add(p);rows.append(f'{sha256_file(p)}  {p}\n')
    (out/'STAGE05F2E_SHA256.tsv').write_text(''.join(rows))
    print(json.dumps(gate,indent=2))

def selftest():
    # Translation and coordinates.
    s='ATGAAA'*25+'TAA'+'ATGCCC'*25+'TAG'
    seg=list(iter_whole_genome_segments('X',s))
    assert any(x['strand']=='+' and x['aa_length']>=45 for x in seg)
    # thresholds
    r={'evalue':'1e-8','qcov':'0.9','fident':'0.3'}
    assert hit_class(r,False)=='CORE_SEQUENCE_HOMOLOG'
    assert longer(100,130) and not major(100,130) and major(100,180)
    # union-find
    uf=UF();uf.union('a','b');assert uf.find('a')==uf.find('b')
    print('[PASS] Stage05F2E translation/threshold/component selftest')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--selftest',action='store_true')
    ap.add_argument('--f2a');ap.add_argument('--f2d');ap.add_argument('--project-root')
    ap.add_argument('--threads',type=int,default=32)
    ap.add_argument('--reuse-v100-work',default='')
    ap.add_argument('--outdir',default='outputs');ap.add_argument('--rawdir',default='raw');ap.add_argument('--work',default='work');ap.add_argument('--logs',default='logs')
    a=ap.parse_args()
    if a.selftest:selftest();return
    if not (a.f2a and a.f2d and a.project_root):ap.error('--f2a --f2d --project-root required unless --selftest')
    run(a)
if __name__=='__main__':main()
