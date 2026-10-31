#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, re, tarfile
from collections import Counter, defaultdict
from pathlib import Path

EXPECTED_QUERY_N=30
EXPECTED_F1_UNIQUE_N=380969
EXPECTED_F1_INSTANCE_N=498879
SEARCH_E=1e-2
CORE_E=1e-5
CORE_QCOV=0.70
CORE_FIDENT=0.25
BROAD_E=1e-3
BROAD_QCOV=0.50
MAX_SEQS=10000

F1_NAMES=[
 'STAGE05F1_UNIQUE_SIXFRAME_PEPTIDES.faa',
 'STAGE05F1_SIXFRAME_PEPTIDE_INSTANCES.tsv',
 'STAGE05F1_RESCUE_PRIORITY.tsv',
 'STAGE05F1_GATE.json']
F2A_NAMES=[
 'STAGE05F2A_NOVELTY_FIRST_STRICT.faa',
 'STAGE05F2A_NOVELTY_FIRST_STRICT_RESCUES.tsv',
 'STAGE05F2A_GATE.json']
F2CB_NAMES=['STAGE05F2CB_SUMMARY.json']

def read_tsv(p):
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows,cols=None):
    p.parent.mkdir(parents=True,exist_ok=True)
    if cols is None: cols=list(rows[0]) if rows else []
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore')
        if cols: w.writeheader()
        w.writerows(rows)

def asbool(x): return str(x).strip().lower() in ('1','true','yes','y')
def ff(x,default=float('nan')):
    try:return float(x)
    except:return default

def splitsemi(x): return [z for z in str(x or '').split(';') if z]

def sha(s): return hashlib.sha256(s.encode()).hexdigest()

def resolve_dir_or_tar(x:Path,names:list[str],work:Path,tag:str):
    if x.is_dir():
        out=x/'outputs'
        if all((out/n).exists() for n in names): return out,'live_directory'
    if x.is_file() and x.name.endswith(('.tar.gz','.tgz')):
        d=work/tag; d.mkdir(parents=True,exist_ok=True)
        wanted=set(names)
        with tarfile.open(x,'r:gz') as tf:
            for m in tf.getmembers():
                if m.isfile() and Path(m.name).name in wanted:
                    fh=tf.extractfile(m)
                    if fh: (d/Path(m.name).name).write_bytes(fh.read())
        if all((d/n).exists() for n in names): return d,'result_tar'
    raise RuntimeError(f'cannot resolve {names} from {x}')

def parse_fasta(path:Path, keep_ids=None):
    # returns dict id -> {seq,sha,len,header}; if keep_ids None, all records
    out={}; h=None; seq=[]; raw=''
    def emit():
        if h is None:return
        if keep_ids is None or h in keep_ids:
            s=''.join(seq).replace('*','').upper()
            out[h]={'seq':s,'sha256':sha(s),'aa_length':len(s),'header':raw}
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line:continue
            if line.startswith('>'):
                emit(); raw=line[1:]; h=raw.split()[0]; seq=[]
            else: seq.append(line)
        emit()
    return out

def count_fasta(path:Path):
    n=0
    with open(path) as f:
        for line in f:
            if line.startswith('>'): n+=1
    return n

def count_tsv_data(path:Path):
    with open(path) as f:
        n=-1
        for n,_ in enumerate(f): pass
    return max(0,n)

def parse_search(path:Path):
    fields=['query','target','evalue','bits','fident','alnlen','qcov','tcov','qlen','tlen']
    if not path.exists() or path.stat().st_size==0:return []
    with open(path) as f:
        return list(csv.DictReader(f,delimiter='\t',fieldnames=fields))

def hit_class(r, exact):
    if exact:return 'EXACT_AA'
    e=ff(r['evalue'],9e99); qc=ff(r['qcov'],0); fi=ff(r['fident'],0)
    if e<=CORE_E and qc>=CORE_QCOV and fi>=CORE_FIDENT:return 'CORE_SEQUENCE_HOMOLOG'
    if e<=BROAD_E and qc>=BROAD_QCOV:return 'BROAD_SEQUENCE_HOMOLOG_SENSITIVITY'
    return 'SEARCH_RESERVOIR'

def longer(qlen,tlen): return tlen>=qlen+30 and tlen>=1.25*qlen
def major_ext(qlen,tlen): return longer(qlen,tlen) and tlen>=max(180,1.5*qlen)

def prepare(args):
    work=Path(args.work); work.mkdir(parents=True,exist_ok=True)
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    f1,r1=resolve_dir_or_tar(Path(args.f1).resolve(),F1_NAMES,work,'f1')
    f2a,r2=resolve_dir_or_tar(Path(args.f2a).resolve(),F2A_NAMES,work,'f2a')
    g1=json.loads((f1/'STAGE05F1_GATE.json').read_text()); g2=json.loads((f2a/'STAGE05F2A_GATE.json').read_text())
    if not g1.get('gate_pass') or not g2.get('gate_pass'): raise RuntimeError('upstream gate not PASS')
    qrows=read_tsv(f2a/'STAGE05F2A_NOVELTY_FIRST_STRICT_RESCUES.tsv')
    if len(qrows)!=EXPECTED_QUERY_N: raise RuntimeError(f'expected {EXPECTED_QUERY_N} novelty queries, got {len(qrows)}')
    qfa=f2a/'STAGE05F2A_NOVELTY_FIRST_STRICT.faa'; qseq=parse_fasta(qfa)
    if len(qseq)!=EXPECTED_QUERY_N: raise RuntimeError(f'query FASTA count {len(qseq)} != {EXPECTED_QUERY_N}')
    q_by_id={r['query_id']:r for r in qrows}
    if set(qseq)!=set(q_by_id): raise RuntimeError('query IDs differ between F2A TSV and FASTA')
    for q,x in qseq.items():
        if x['sha256']!=q_by_id[q]['aa_sha256']: raise RuntimeError(f'F2A hash mismatch {q}')
    tf=f1/'STAGE05F1_UNIQUE_SIXFRAME_PEPTIDES.faa'; inst=f1/'STAGE05F1_SIXFRAME_PEPTIDE_INSTANCES.tsv'
    tn=count_fasta(tf); inn=count_tsv_data(inst)
    if tn!=EXPECTED_F1_UNIQUE_N: raise RuntimeError(f'expected F1 unique {EXPECTED_F1_UNIQUE_N}, got {tn}')
    if inn!=EXPECTED_F1_INSTANCE_N: raise RuntimeError(f'expected F1 instances {EXPECTED_F1_INSTANCE_N}, got {inn}')
    # exact query IDs/hashes must be present in target F1 universe
    f1q=parse_fasta(tf,keep_ids=set(qseq))
    if set(f1q)!=set(qseq): raise RuntimeError(f'{len(set(qseq)-set(f1q))} F2A queries absent from F1 unique FASTA')
    for q in qseq:
        if qseq[q]['sha256']!=f1q[q]['sha256']: raise RuntimeError(f'F1/F2A query hash mismatch {q}')
    f2cb_status='NOT_SUPPLIED'; f2cb_route='none'; f2cb_query=''; f2cb_tier=''
    if args.f2cb and args.f2cb.upper()!='NONE':
        c,rc=resolve_dir_or_tar(Path(args.f2cb).resolve(),F2CB_NAMES,work,'f2cb')
        s=json.loads((c/'STAGE05F2CB_SUMMARY.json').read_text())
        if not s.get('gate_pass'): raise RuntimeError('F2CB gate not PASS')
        f2cb_route=rc; f2cb_query=str(s.get('query','')); f2cb_tier=str(s.get('best_tier',''))
        f2cb_status=f'{f2cb_query}:{f2cb_tier}'
    prep={
      'stage':'Stage05F2D','version':'v1.0.0','phase':'PREPARE','gate':'PASS_STAGE05F2D_FROZEN_FAMILY_EXPANSION_PREPARED','gate_pass':True,
      'f1_route':r1,'f2a_route':r2,'f2cb_route':f2cb_route,'novelty_query_n':len(qrows),'f1_unique_target_n':tn,'f1_instance_n':inn,
      'f2cb_control_status':f2cb_status,'search_evalue':SEARCH_E,'core_homolog_rule':f'e<={CORE_E};qcov>={CORE_QCOV};fident>={CORE_FIDENT}',
      'broad_rule':f'e<={BROAD_E};qcov>={BROAD_QCOV}','max_seqs_per_query':MAX_SEQS,
      'candidate_membership_changed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False}
    (out/'STAGE05F2D_PREPARE_GATE.json').write_text(json.dumps(prep,indent=2)+'\n')
    # Copy exact frozen query FASTA into work for search provenance.
    (work/'STAGE05F2D_QUERY_30.faa').write_bytes(qfa.read_bytes())
    paths={'f1_outputs':str(f1),'f2a_outputs':str(f2a),'target_fasta':str(tf),'instance_tsv':str(inst),'query_fasta':str(work/'STAGE05F2D_QUERY_30.faa')}
    (work/'RESOLVED_INPUTS.json').write_text(json.dumps(paths,indent=2)+'\n')
    print(json.dumps(prep,indent=2))

def summarize(args):
    out=Path(args.outdir); raw=Path(args.m8); work=Path(args.work)
    paths=json.loads((work/'RESOLVED_INPUTS.json').read_text())
    f1=Path(paths['f1_outputs']); f2a=Path(paths['f2a_outputs']); tf=Path(paths['target_fasta']); instp=Path(paths['instance_tsv'])
    qrows=read_tsv(f2a/'STAGE05F2A_NOVELTY_FIRST_STRICT_RESCUES.tsv'); qmeta={r['query_id']:r for r in qrows}
    qseq=parse_fasta(Path(paths['query_fasta']))
    hits=parse_search(raw)
    byq=defaultdict(list)
    target_ids=set()
    for r in hits:
        byq[r['query']].append(r); target_ids.add(r['target'])
    # Resolve only targets actually returned by MMseqs.
    tseq=parse_fasta(tf,keep_ids=target_ids)
    missing=target_ids-set(tseq)
    if missing: raise RuntimeError(f'missing {len(missing)} search targets in F1 FASTA, first={sorted(missing)[:10]}')
    # Carry forward F1 HMM/rescue metadata for returned targets.
    tprior={}
    wanted=target_ids
    with open(f1/'STAGE05F1_RESCUE_PRIORITY.tsv',newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r.get('query_id') in wanted:tprior[r['query_id']]=r
    # Annotate search rows and identify relevant hashes.
    annot=[]; relevant_hashes=set(); longer_hashes=set(); target_to_queries=defaultdict(set)
    exact_recovered=set(); cap=Counter()
    for q in qseq:
        cap[q]=len(byq.get(q,[]))
        for r in byq.get(q,[]):
            t=r['target']; tx=tseq[t]; exact=(tx['sha256']==qseq[q]['sha256'])
            hc=hit_class(r,exact); ql=int(float(r['qlen'])); tl=int(float(r['tlen']))
            lg=(not exact and hc=='CORE_SEQUENCE_HOMOLOG' and longer(ql,tl)); mj=(lg and major_ext(ql,tl))
            if exact: exact_recovered.add(q)
            if hc in ('EXACT_AA','CORE_SEQUENCE_HOMOLOG','BROAD_SEQUENCE_HOMOLOG_SENSITIVITY'):
                relevant_hashes.add(tx['sha256']); target_to_queries[tx['sha256']].add(q)
            if lg: longer_hashes.add(tx['sha256'])
            pr=tprior.get(t,{})
            x=dict(r); x.update({
              'query_aa_sha256':qseq[q]['sha256'],'target_aa_sha256':tx['sha256'],'hit_class':hc,'exact_aa':str(exact),
              'longer_core_homolog':str(lg),'major_extension':str(mj),
              'target_formal_rvt_profile':pr.get('formal_rvt_profile',''),'target_exploratory_rvt_profile':pr.get('exploratory_rvt_profile',''),
              'target_best_rvt_target':pr.get('best_rvt_target',''),'target_best_full_evalue':pr.get('best_full_evalue',''),
              'target_rescue_class':pr.get('rescue_class',''),'target_prodigal_relations':pr.get('prodigal_relations','')})
            annot.append(x)
    if len(exact_recovered)!=EXPECTED_QUERY_N: raise RuntimeError(f'exact query self-recovery {len(exact_recovered)}/{EXPECTED_QUERY_N}')
    # Stream instance registry only once; aggregate accessions/candidates/relations by hash.
    agg={h:{'accessions':set(),'candidate_keys':set(),'relations':set(),'instance_n':0,'edge_partial_n':0} for h in relevant_hashes}
    longer_instances=[]
    with open(instp,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            h=r.get('aa_sha256','')
            if h not in agg: continue
            a=agg[h]; a['instance_n']+=1
            if r.get('sequence_id'):a['accessions'].add(r['sequence_id'])
            a['candidate_keys'].update(splitsemi(r.get('candidate_keys','')))
            if r.get('prodigal_relation'):a['relations'].add(r['prodigal_relation'])
            if asbool(r.get('edge_partial')):a['edge_partial_n']+=1
            if h in longer_hashes: longer_instances.append(r)
    # Enrich hit table with instance counts.
    for x in annot:
        a=agg.get(x['target_aa_sha256'],{'accessions':set(),'candidate_keys':set(),'relations':set(),'instance_n':0,'edge_partial_n':0})
        x.update({'target_instance_n':a['instance_n'],'target_genome_n':len(a['accessions']),'target_accessions':';'.join(sorted(a['accessions'])),
                  'target_candidate_n':len(a['candidate_keys']),'target_candidate_keys':';'.join(sorted(a['candidate_keys'])),
                  'target_instance_relations':';'.join(sorted(a['relations'])),'target_edge_partial_instance_n':a['edge_partial_n']})
    cols=list(annot[0]) if annot else []
    write_tsv(out/'STAGE05F2D_ALL_SEARCH_HITS_ANNOTATED.tsv',annot,cols)
    # Query-level summaries.
    summaries=[]; ext_rows=[]; priority_rows=[]; major_target_ids=set()
    priority_counts=Counter()
    for q in sorted(qseq):
        rs=[x for x in annot if x['query']==q]
        original=set(splitsemi(qmeta[q].get('accessions','')))
        exact=[x for x in rs if x['hit_class']=='EXACT_AA']
        core=[x for x in rs if x['hit_class']=='CORE_SEQUENCE_HOMOLOG']
        broad=[x for x in rs if x['hit_class']=='BROAD_SEQUENCE_HOMOLOG_SENSITIVITY']
        exts=[x for x in core if asbool(x['longer_core_homolog'])]
        majors=[x for x in exts if asbool(x['major_extension'])]
        fam_acc=set(); core_acc=set(); broad_acc=set(); exact_acc=set(); major_acc=set()
        for z in exact: exact_acc.update(splitsemi(z.get('target_accessions','')))
        for z in core: core_acc.update(splitsemi(z.get('target_accessions','')))
        for z in broad: broad_acc.update(splitsemi(z.get('target_accessions','')))
        for z in majors: major_acc.update(splitsemi(z.get('target_accessions','')))
        fam_acc=exact_acc|core_acc|broad_acc
        formal_major=[z for z in majors if asbool(z.get('target_formal_rvt_profile'))]
        if formal_major and len(major_acc)>=2: pri='P1_LONGER_FORMAL_RVT_RESCUE_MULTI_GENOME'
        elif formal_major: pri='P2_LONGER_FORMAL_RVT_RESCUE_SINGLE_GENOME'
        elif majors: pri='P3_MAJOR_EXTENSION_NO_FORMAL_RVT'
        elif len(core_acc)>=2: pri='P4_CORE_REPLICATED_NO_MAJOR_EXTENSION'
        else: pri='P5_NO_CORE_FAMILY_EXPANSION'
        priority_counts[pri]+=1
        best=min(core,key=lambda z:ff(z['evalue'],9e99)) if core else None
        longest=max(core,key=lambda z:int(float(z['tlen']))) if core else None
        control='FROZEN_RT_FOLD_NO_HIT_CONTROL' if q=='SF_00172800' else ''
        s={
          'query_id':q,'query_aa_sha256':qseq[q]['sha256'],'query_aa_length':qseq[q]['aa_length'],'original_genome_n':len(original),'original_accessions':';'.join(sorted(original)),
          'search_hit_n':len(rs),'search_cap_reached':str(len(rs)>=MAX_SEQS),'exact_replication_genome_n':len(exact_acc),'exact_replication_accessions':';'.join(sorted(exact_acc)),
          'core_homolog_unique_aa_n':len({z['target_aa_sha256'] for z in core}),'core_homolog_genome_n':len(core_acc),'core_new_genome_n':len(core_acc-original),
          'broad_homolog_unique_aa_n':len({z['target_aa_sha256'] for z in broad}),'broad_homolog_genome_n':len(broad_acc),
          'family_total_genome_n':len(fam_acc),'longer_core_unique_aa_n':len({z['target_aa_sha256'] for z in exts}),'major_extension_unique_aa_n':len({z['target_aa_sha256'] for z in majors}),'major_extension_genome_n':len(major_acc),
          'formal_rvt_major_extension_unique_aa_n':len({z['target_aa_sha256'] for z in formal_major}),'best_core_target':best['target'] if best else '',
          'best_core_evalue':best['evalue'] if best else '','best_core_qcov':best['qcov'] if best else '','best_core_fident':best['fident'] if best else '',
          'longest_core_target':longest['target'] if longest else '','longest_core_target_len':longest['tlen'] if longest else '',
          'max_core_length_ratio':(max([int(float(z['tlen']))/qseq[q]['aa_length'] for z in core]) if core else ''),
          'f2d_priority':pri,'f2cb_control_status':control}
        summaries.append(s)
        for z in exts:
            y=dict(z); y['f2d_priority']=pri; ext_rows.append(y)
            if asbool(z['major_extension']): major_target_ids.add(z['target'])
        if pri!='P5_NO_CORE_FAMILY_EXPANSION': priority_rows.append(s)
    write_tsv(out/'STAGE05F2D_FAMILY_SUMMARY.tsv',summaries)
    write_tsv(out/'STAGE05F2D_PRIORITY_FAMILIES.tsv',priority_rows)
    write_tsv(out/'STAGE05F2D_LONGER_ORF_CANDIDATES.tsv',ext_rows, list(ext_rows[0]) if ext_rows else cols+['f2d_priority'])
    write_tsv(out/'STAGE05F2D_LONGER_ORF_INSTANCE_MAP.tsv',longer_instances)
    # FASTA for all major-extension unique targets.
    with (out/'STAGE05F2D_MAJOR_EXTENSION_TARGETS.faa').open('w') as f:
        for t in sorted(major_target_ids):
            x=tseq[t]
            pr=tprior.get(t,{})
            f.write(f'>{t} sha256={x["sha256"]} len={x["aa_length"]} formal_rvt={pr.get("formal_rvt_profile","")} profile={pr.get("best_rvt_target","")} rescue={pr.get("rescue_class","")}\n{x["seq"]}\n')
    prep=json.loads((out/'STAGE05F2D_PREPARE_GATE.json').read_text())
    gate={
      'stage':'Stage05F2D','version':'v1.0.0','gate':'PASS_STAGE05F2D_ANNOTATION_INDEPENDENT_FAMILY_EXPANSION','gate_pass':True,
      'f1_route':prep['f1_route'],'f2a_route':prep['f2a_route'],'f2cb_control_status':prep['f2cb_control_status'],
      'novelty_query_n':EXPECTED_QUERY_N,'f1_unique_target_n':EXPECTED_F1_UNIQUE_N,'f1_instance_n':EXPECTED_F1_INSTANCE_N,
      'mmseqs_search_hit_rows_n':len(annot),'exact_query_recovery_n':len(exact_recovered),'search_cap_reached_query_n':sum(len(byq.get(q,[]))>=MAX_SEQS for q in qseq),
      'core_expansion_query_n':sum(int(s['core_homolog_unique_aa_n'])>0 for s in summaries),'major_extension_query_n':sum(int(s['major_extension_unique_aa_n'])>0 for s in summaries),
      'formal_rvt_major_extension_query_n':sum(int(s['formal_rvt_major_extension_unique_aa_n'])>0 for s in summaries),'priority_counts':dict(priority_counts),
      'new_hmm_search_performed':False,'new_structure_prediction_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'candidate_membership_changed':False,
      'claim_boundary':'Family expansion and longer annotation-independent stop-to-stop peptide recovery only. Homology, replication, inherited F1 profile status, and extension do not establish biochemical RT activity or ART identity.'}
    (out/'STAGE05F2D_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    notes=[
      '# Stage05F2D manuscript-facing notes\n\n',
      f'- Frozen novelty-first query sequences: {EXPECTED_QUERY_N}.\n',
      f'- Stage05F1 annotation-independent target universe: {EXPECTED_F1_UNIQUE_N:,} unique AA / {EXPECTED_F1_INSTANCE_N:,} genomic instances.\n',
      f'- Queries with core non-exact family expansion: {gate["core_expansion_query_n"]}.\n',
      f'- Queries with major longer-peptide extension: {gate["major_extension_query_n"]}.\n',
      f'- Queries whose major extension already carried a formal F1 RVT profile: {gate["formal_rvt_major_extension_query_n"]}.\n\n',
      'The search universe remains annotation-independent at peptide generation. Existing Prodigal relations and RVT HMM calls are carried forward only for interpretation. SF_00172800, if present, remains a frozen RT-fold NO_HIT control rather than a biological negative.\n']
    (out/'STAGE05F2D_MANUSCRIPT_NOTES.md').write_text(''.join(notes))
    print(json.dumps(gate,indent=2))

def main():
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='cmd',required=True)
    p=sub.add_parser('prepare'); p.add_argument('--f1',required=True);p.add_argument('--f2a',required=True);p.add_argument('--f2cb',default='NONE');p.add_argument('--outdir',default='outputs');p.add_argument('--work',default='work')
    s=sub.add_parser('summarize');s.add_argument('--m8',required=True);s.add_argument('--outdir',default='outputs');s.add_argument('--work',default='work')
    a=ap.parse_args(); prepare(a) if a.cmd=='prepare' else summarize(a)
if __name__=='__main__':main()
