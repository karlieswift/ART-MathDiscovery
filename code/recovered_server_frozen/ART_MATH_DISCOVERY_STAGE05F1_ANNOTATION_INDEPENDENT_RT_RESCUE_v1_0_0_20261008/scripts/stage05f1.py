#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, os, re, shutil, subprocess, sys, tarfile
from collections import defaultdict
from pathlib import Path

EXPECTED_DUAL_SHA = "04706102bdb8129e105c6cdba915d33e95e19963463869552622ad73283729ae"
EXPECTED_SOURCE_MANIFEST_SHA = "5eec6d7ee00584c6266b8a953b6112540148b949341e400099b285ce30d235f0"
FORMAL_E = 1e-4
EXPLORE_E = 1e-2
MIN_AA = 45
FLANK = 10000

CODON_TABLE_11 = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
RC = str.maketrans('ACGTN','TGCAN')

def sha256_file(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def open_text(p: Path):
    return gzip.open(p,'rt') if p.suffix=='.gz' else p.open('r', newline='')

def read_fasta_one(p: Path):
    name=None; seq=[]
    with open_text(p) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: break
                name=line[1:].split()[0]
            else: seq.append(line.upper())
    if name is None: raise RuntimeError(f'no FASTA record in {p}')
    return name,''.join(seq)

def read_tsv(p: Path):
    with open_text(p) as f:
        return list(csv.DictReader(f, delimiter='\t'))

def header(p: Path):
    with open_text(p) as f:
        r=csv.reader(f, delimiter='\t')
        return next(r, [])

def candidate_files(root: Path, names):
    for pat in names:
        for p in root.rglob(pat):
            if p.is_file(): yield p

def resolve_by_sha(root: Path, expected: str, patterns):
    override=os.environ.get('STAGE05F1_DUAL_TSV' if expected==EXPECTED_DUAL_SHA else 'STAGE05F1_SOURCE_MANIFEST','').strip()
    if override:
        p=Path(override)
        if sha256_file(p)!=expected: raise RuntimeError(f'override hash mismatch {p}')
        return p,'override'
    seen=set()
    for p in candidate_files(root, patterns):
        rp=str(p.resolve())
        if rp in seen: continue
        seen.add(rp)
        try:
            if p.stat().st_size > 200_000_000: continue
            if sha256_file(p)==expected: return p,'live_exact_sha256'
        except Exception: pass
    # tar fallback: scan only Stage05B/C1 archives, extract exact byte match
    for p in sorted(root.glob('**/*STAGE05B*.tar.gz')) + sorted(root.glob('**/*STAGE05C1*.tar.gz')):
        try:
            with tarfile.open(p,'r:gz') as tf:
                for m in tf.getmembers():
                    if not m.isfile() or m.size>200_000_000: continue
                    bn=Path(m.name).name.lower()
                    if not any(x in bn for x in ('candidate','dual','manifest','region')): continue
                    fh=tf.extractfile(m)
                    if fh is None: continue
                    data=fh.read(); h=hashlib.sha256(data).hexdigest()
                    if h==expected:
                        out=root/'.stage05f1_resolved'/Path(m.name).name
                        out.parent.mkdir(parents=True,exist_ok=True); out.write_bytes(data)
                        return out,f'tar:{p.name}:{m.name}'
        except Exception: pass
    raise RuntimeError(f'could not resolve frozen file sha256={expected}')

def detect_cols(rows, kind='candidate'):
    if not rows: raise RuntimeError('empty table')
    cols=list(rows[0])
    def choose(opts, required=True):
        for x in opts:
            if x in cols: return x
        if required: raise RuntimeError(f'missing columns {opts}; have {cols[:50]}')
        return None
    acc=choose(['sequence_id','accession','seq_id','genome_id'])
    if kind=='candidate':
        st=choose(['start0','region_start0','union_start0','candidate_start0','start'])
        en=choose(['end0','region_end0','union_end0','candidate_end0','end'])
        rid=choose(['region_id','candidate_id','union_region_id','id'],False)
        rank=choose(['rank','region_rank','priority_rank'],False)
        return acc,st,en,rid,rank
    else:
        st=choose(['start0','cds_start0','orf_start0','protein_start0','start'])
        en=choose(['end0','cds_end0','orf_end0','protein_end0','end'])
        strand=choose(['strand','cds_strand','orf_strand'],False)
        key=choose(['record_key','protein_key','orf_id','target_name','id'],False)
        return acc,st,en,strand,key

def resolve_fasta_manifest(root: Path, manifest: Path):
    rows=read_tsv(manifest); cols=list(rows[0]) if rows else []
    acc_col=next((x for x in ['sequence_id','accession','seq_id','genome_id'] if x in cols),None)
    path_col=next((x for x in cols if 'fasta' in x.lower() and ('path' in x.lower() or 'file' in x.lower() or 'frozen' in x.lower())),None)
    if not acc_col or not path_col: raise RuntimeError(f'cannot infer fasta columns from {cols}')
    out={}
    for r in rows:
        acc=r[acc_col].strip(); q=r[path_col].strip()
        if not acc or not q: continue
        cands=[Path(q), manifest.parent/q, root/q]
        # basename fallback
        cands += list(root.glob(f'**/{Path(q).name}'))[:3]
        hit=next((p for p in cands if p.exists() and p.is_file()),None)
        if hit: out[acc]=hit.resolve()
    return out,acc_col,path_col,len(rows)

def resolve_hmm(root: Path):
    override=os.environ.get('STAGE05F1_RVT_HMM','').strip()
    if override and Path(override).exists(): return Path(override),'override'
    cands=[]
    for base in root.glob('ART_MATH_DISCOVERY_STAGE05D1*'):
        for p in base.rglob('*.hmm'):
            try:
                txt=p.read_text(errors='ignore')
                n=txt.count('HMMER3/f')
                if 40<=n<=60: cands.append((abs(n-45),p,n))
            except Exception: pass
    if not cands:
        for p in root.rglob('*.hmm'):
            if 'RVT' not in p.name.upper() and 'RT' not in p.name.upper(): continue
            try:
                txt=p.read_text(errors='ignore'); n=txt.count('HMMER3/f')
                if 40<=n<=60: cands.append((abs(n-45),p,n))
            except Exception: pass
    if not cands: raise RuntimeError('could not resolve local ~45-profile RVT-All HMM; set STAGE05F1_RVT_HMM')
    cands.sort(key=lambda x:(x[0],len(str(x[1]))))
    return cands[0][1],f'auto_models={cands[0][2]}'

def resolve_orf_registry(root: Path):
    override=os.environ.get('STAGE05F1_ORF_REGISTRY','').strip()
    if override and Path(override).exists(): return Path(override),'override'
    bases=list(root.glob('ART_MATH_DISCOVERY_STAGE05C1*'))
    candidates=[]
    for b in bases:
        for p in b.rglob('*.tsv'):
            try:
                h=header(p); low=[x.lower() for x in h]
                if any(x in low for x in ['sequence_id','accession']) and any('start0'==x or 'cds_start0'==x or 'orf_start0'==x for x in low) and any('end0'==x or 'cds_end0'==x or 'orf_end0'==x for x in low):
                    score=(5 if 'protein' in p.name.lower() or 'orf' in p.name.lower() or 'registry' in p.name.lower() else 0) + min(p.stat().st_size/1e6,20)
                    candidates.append((score,p))
            except Exception: pass
    if not candidates: return None,'not_resolved'
    candidates.sort(reverse=True,key=lambda x:x[0]); return candidates[0][1],'auto'

def merge_intervals(items):
    items=sorted(items)
    out=[]
    for s,e in items:
        if not out or s>out[-1][1]: out.append([s,e])
        else: out[-1][1]=max(out[-1][1],e)
    return [(a,b) for a,b in out]

def translate_frame(nt):
    aa=[]
    for i in range(0,len(nt)-2,3): aa.append(CODON_TABLE_11.get(nt[i:i+3],'X'))
    return ''.join(aa)

def segments_from_interval(acc, seq, ws, we):
    sub=seq[ws:we]
    for strand in ['+','-']:
        nt=sub if strand=='+' else sub.translate(RC)[::-1]
        L=len(sub)
        for frame in (0,1,2):
            aa=translate_frame(nt[frame:])
            aa0=0
            for part in aa.split('*'):
                aa1=aa0+len(part)
                if len(part)>=MIN_AA:
                    # nucleotide coords within oriented sub
                    n0=frame+aa0*3; n1=frame+aa1*3
                    if strand=='+': gs,ge=ws+n0,ws+n1
                    else: gs,ge=we-n1,we-n0
                    edge_partial=(aa0==0 or aa1==len(aa))
                    yield {'sequence_id':acc,'start0':gs,'end0':ge,'strand':strand,'frame':frame,'aa_sequence':part,'aa_length':len(part),'edge_partial':edge_partial}
                aa0=aa1+1

def parse_domtbl(p: Path):
    rows=[]
    if not p.exists(): return rows
    with p.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            x=line.rstrip('\n').split(maxsplit=22)
            if len(x)<23: continue
            rows.append({'target_name':x[0],'query_name':x[3],'full_evalue':float(x[6]),'domain_i_evalue':float(x[12]),'domain_score':float(x[13]),'ali_from':int(x[17]),'ali_to':int(x[18]),'description':x[22]})
    return rows

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('project_root'); ap.add_argument('--threads',type=int,default=16); ap.add_argument('--outdir',default='outputs')
    a=ap.parse_args(); root=Path(a.project_root).resolve(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); work=Path('work'); work.mkdir(exist_ok=True)
    dual,dual_route=resolve_by_sha(root,EXPECTED_DUAL_SHA,['*DUAL*.tsv','*dual*.tsv','*CANDIDATE*.tsv','*candidate*.tsv','*REGION*.tsv'])
    mani,mani_route=resolve_by_sha(root,EXPECTED_SOURCE_MANIFEST_SHA,['*MANIFEST*.tsv','*manifest*.tsv'])
    candidates=read_tsv(dual); accc,stc,enc,ridc,rankc=detect_cols(candidates,'candidate')
    if len(candidates)!=3674: raise RuntimeError(f'expected 3674 dual regions, got {len(candidates)}')
    fasta_map,ma,mp,mn=resolve_fasta_manifest(root,mani)
    accessions=sorted(set(r[accc] for r in candidates))
    missing=[x for x in accessions if x not in fasta_map]
    if missing: raise RuntimeError(f'missing FASTA for {len(missing)} candidate genomes, first={missing[:10]}')
    hmm,hmm_route=resolve_hmm(root)
    orf_path,orf_route=resolve_orf_registry(root)
    orf_by_acc=defaultdict(list)
    if orf_path:
        rows=read_tsv(orf_path); oa,os_,oe,ostr,okey=detect_cols(rows,'orf')
        for r in rows:
            try: orf_by_acc[r[oa]].append((int(float(r[os_])),int(float(r[oe])),r.get(ostr,'') if ostr else '',r.get(okey,'') if okey else ''))
            except: pass
        for x in orf_by_acc: orf_by_acc[x].sort()
    # build candidate metadata + merged windows
    cand_meta=[]; windows_by_acc=defaultdict(list)
    for i,r in enumerate(candidates):
        acc=r[accc]; s=int(float(r[stc])); e=int(float(r[enc])); rid=r[ridc] if ridc else f'R{i+1:04d}'; rank=r[rankc] if rankc else ''
        _,seq=read_fasta_one(fasta_map[acc]); L=len(seq)
        ws=max(0,s-FLANK); we=min(L,e+FLANK)
        ck=f'{acc}::{rid}'
        cand_meta.append((ck,acc,s,e,ws,we,rank,L))
        windows_by_acc[acc].append((ws,we))
    # translate merged windows
    peptide_instances=[]; pep_id=0
    for acc in accessions:
        _,seq=read_fasta_one(fasta_map[acc])
        for ws,we in merge_intervals(windows_by_acc[acc]):
            for p in segments_from_interval(acc,seq,ws,we):
                pep_id+=1; p['peptide_instance_id']=f'R6F_{pep_id:08d}'; p['merged_window_start0']=ws; p['merged_window_end0']=we
                aa=p['aa_sequence']; p['aa_sha256']=hashlib.sha256(aa.encode()).hexdigest(); peptide_instances.append(p)
    # map peptide to candidate neighborhoods
    cands_by_acc=defaultdict(list)
    for row in cand_meta: cands_by_acc[row[1]].append(row)
    p2c=defaultdict(list)
    for p in peptide_instances:
        for ck,acc,s,e,ws,we,rank,L in cands_by_acc[p['sequence_id']]:
            if p['end0']>ws and p['start0']<we: p2c[p['peptide_instance_id']].append(ck)
    # annotation comparison post rescue
    for p in peptide_instances:
        overlaps=[]; same=[]
        for s,e,strand,key in orf_by_acc.get(p['sequence_id'],[]):
            if e<=p['start0']: continue
            if s>=p['end0']: break
            if e>p['start0'] and s<p['end0']:
                overlaps.append((s,e,strand,key))
                if not strand or strand==p['strand']: same.append((s,e,strand,key))
        if not orf_path: cls='ORF_REGISTRY_UNRESOLVED'
        elif not same: cls='NO_SAME_STRAND_PRODIGAL_OVERLAP'
        elif any(s<=p['start0'] and e>=p['end0'] for s,e,_,_ in same): cls='CONTAINED_IN_EXISTING_PRODIGAL_ORF'
        else: cls='BOUNDARY_OR_PARTIAL_DISCORDANT_WITH_PRODIGAL'
        p['prodigal_relation']=cls; p['same_strand_prodigal_overlap_n']=len(same)
    # dedupe by aa hash
    unique={}
    for p in peptide_instances:
        h=p['aa_sha256']
        if h not in unique: unique[h]={'aa_sha256':h,'aa_sequence':p['aa_sequence'],'aa_length':p['aa_length'],'instance_ids':[],'accessions':set(),'candidate_keys':set(),'relations':set(),'motif_yxdd_like':bool(re.search(r'[YF][A-Z]DD',p['aa_sequence']))}
        u=unique[h]; u['instance_ids'].append(p['peptide_instance_id']); u['accessions'].add(p['sequence_id']); u['relations'].add(p['prodigal_relation']);
        for ck in p2c[p['peptide_instance_id']]: u['candidate_keys'].add(ck)
    # FASTA + registries
    fa=out/'STAGE05F1_UNIQUE_SIXFRAME_PEPTIDES.faa'
    id_by_hash={}
    with fa.open('w') as f:
        for j,(h,u) in enumerate(sorted(unique.items()),1):
            q=f'SF_{j:08d}'; id_by_hash[h]=q; f.write(f'>{q} sha256={h} aa={u["aa_length"]} instances={len(u["instance_ids"])} genomes={len(u["accessions"])} candidates={len(u["candidate_keys"])}\n{u["aa_sequence"]}\n')
    with (out/'STAGE05F1_SIXFRAME_PEPTIDE_INSTANCES.tsv').open('w',newline='') as f:
        cols=['peptide_instance_id','sequence_id','start0','end0','strand','frame','aa_length','aa_sha256','edge_partial','prodigal_relation','same_strand_prodigal_overlap_n','candidate_keys']
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t'); w.writeheader()
        for p in peptide_instances:
            q={k:p.get(k,'') for k in cols}; q['candidate_keys']=';'.join(sorted(p2c[p['peptide_instance_id']])) ; w.writerow(q)
    # HMM
    dom=out/'STAGE05F1_RVT_ALL.domtblout'; tbl=out/'STAGE05F1_RVT_ALL.tblout'
    cmd=['hmmscan','--cpu',str(a.threads),'-E',str(EXPLORE_E),'--domE',str(EXPLORE_E),'--domtblout',str(dom),'--tblout',str(tbl),str(hmm),str(fa)]
    subprocess.run(cmd,check=True,stdout=(work/'hmmscan.stdout').open('w'),stderr=(work/'hmmscan.stderr').open('w'))
    hits=parse_domtbl(dom); hits_by_q=defaultdict(list)
    for h in hits: hits_by_q[h['query_name']].append(h)
    # prioritized unique table
    pri=[]
    for h,u in sorted(unique.items()):
        q=id_by_hash[h]; hh=sorted(hits_by_q.get(q,[]),key=lambda x:(x['full_evalue'],x['domain_i_evalue']))
        best=hh[0] if hh else None
        formal=any(x['full_evalue']<=FORMAL_E and x['domain_i_evalue']<=FORMAL_E for x in hh)
        exploratory=bool(hh)
        rels=u['relations']
        gene_caller_rescue=('NO_SAME_STRAND_PRODIGAL_OVERLAP' in rels or 'BOUNDARY_OR_PARTIAL_DISCORDANT_WITH_PRODIGAL' in rels)
        if formal and gene_caller_rescue: cls='FORMAL_PROFILE_GENE_CALLER_RESCUE'
        elif formal: cls='FORMAL_PROFILE_EXISTING_ORF_EXPLAINED_OR_MIXED'
        elif exploratory and gene_caller_rescue: cls='EXPLORATORY_PROFILE_GENE_CALLER_RESCUE'
        elif u['motif_yxdd_like'] and len(u['accessions'])>=2 and gene_caller_rescue: cls='HMM_UNDER_DETECTED_MOTIF_REPLICATED_RESCUE_RESERVOIR'
        elif u['motif_yxdd_like'] and gene_caller_rescue: cls='HMM_UNDER_DETECTED_MOTIF_SINGLETON_RESCUE_RESERVOIR'
        else: cls='BACKGROUND_RESCUED_PEPTIDE'
        pri.append({'query_id':q,'aa_sha256':h,'aa_length':u['aa_length'],'genome_n':len(u['accessions']),'candidate_n':len(u['candidate_keys']),'instance_n':len(u['instance_ids']),'prodigal_relations':';'.join(sorted(rels)),'yxdd_like':u['motif_yxdd_like'],'formal_rvt_profile':formal,'exploratory_rvt_profile':exploratory,'best_rvt_target':best['target_name'] if best else '','best_full_evalue':best['full_evalue'] if best else '','best_domain_i_evalue':best['domain_i_evalue'] if best else '','rescue_class':cls,'candidate_keys':';'.join(sorted(u['candidate_keys'])),'accessions':';'.join(sorted(u['accessions']))})
    with (out/'STAGE05F1_RESCUE_PRIORITY.tsv').open('w',newline='') as f:
        cols=list(pri[0]) if pri else [] ; w=csv.DictWriter(f,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows(pri)
    # optional MMseqs on non-background priority
    priority_ids={r['query_id'] for r in pri if r['rescue_class']!='BACKGROUND_RESCUED_PEPTIDE'}
    pfa=work/'priority.faa'
    seq_by_id={id_by_hash[h]:u['aa_sequence'] for h,u in unique.items()}
    with pfa.open('w') as f:
        for q in sorted(priority_ids): f.write(f'>{q}\n{seq_by_id[q]}\n')
    mmseqs_status='NOT_RUN_NO_PRIORITY_OR_TOOL'
    if priority_ids and shutil.which('mmseqs'):
        try:
            subprocess.run(['mmseqs','easy-cluster',str(pfa),str(work/'mmseqs_cluster'),str(work/'mmseqs_tmp'),'--min-seq-id','0.3','-c','0.5','--cov-mode','0','--threads',str(a.threads)],check=True,stdout=(work/'mmseqs.stdout').open('w'),stderr=(work/'mmseqs.stderr').open('w'))
            src=work/'mmseqs_cluster_cluster.tsv'
            if src.exists(): shutil.copy2(src,out/'STAGE05F1_PRIORITY_MMSEQS_CLUSTER.tsv')
            mmseqs_status='PASS'
        except Exception as e: mmseqs_status=f'FAILED_NONFATAL:{e}'
    # summary
    class_counts=defaultdict(int)
    for r in pri: class_counts[r['rescue_class']]+=1
    summary={'stage':'Stage05F1','version':'v1.0.0','gate':'PASS_STAGE05F1_ANNOTATION_INDEPENDENT_SIXFRAME_RT_RESCUE','gate_pass':True,'dual_candidate_sha256':sha256_file(dual),'dual_candidate_route':dual_route,'source_manifest_sha256':sha256_file(mani),'source_manifest_route':mani_route,'candidate_n':len(candidates),'candidate_genome_n':len(accessions),'flank_bp':FLANK,'translation_mode':'six-frame maximal stop-to-stop; no start codon required; table11; min45aa','sixframe_peptide_instance_n':len(peptide_instances),'sixframe_unique_aa_n':len(unique),'orf_registry':str(orf_path) if orf_path else None,'orf_registry_route':orf_route,'rvt_hmm':str(hmm),'rvt_hmm_route':hmm_route,'formal_hmm_threshold':FORMAL_E,'exploratory_hmm_threshold':EXPLORE_E,'rescue_class_counts':dict(class_counts),'mmseqs_priority_cluster_status':mmseqs_status,'dna_score_recomputed':False,'dna_rank_recomputed':False,'candidate_membership_changed':False,'annotation_used_to_generate_rescue_peptides':False,'claim_boundary':'Stage05F1 identifies annotation-independent peptide rescue evidence only; it does not establish biochemical RT or ART identity.'}
    (out/'STAGE05F1_GATE.json').write_text(json.dumps(summary,indent=2)+'\n')
    # manuscript note
    notes=['# Stage05F1 annotation-independent RT rescue\n','Primary scan: all 3,674 frozen dual-arm regions, ±10 kb, six-frame maximal stop-to-stop translation with no start-codon requirement. Existing gene models were used only after peptide generation to classify whether a segment was already represented.\n','## Rescue class counts\n']
    for k,v in sorted(class_counts.items()): notes.append(f'- {k}: {v}\n')
    notes.append('\nNo HMM or motif hit is sufficient for RT identity. HMM-negative motif-bearing replicated segments remain an unresolved structure/conservation reservoir.\n')
    (out/'STAGE05F1_MANUSCRIPT_NOTES.md').write_text(''.join(notes))
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
