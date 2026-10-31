#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,tarfile,hashlib
from pathlib import Path
from collections import defaultdict,Counter

def read_tsv(p):
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,cols=None):
    if cols is None: cols=list(rows[0]) if rows else []
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
def splitsemi(x): return [z for z in str(x or '').split(';') if z]
def asbool(x): return str(x).lower() in ('true','1','yes')
def flt(x,default=0.0):
    try:return float(x)
    except:return default

def resolve_dir_or_tar(x:Path, names:list[str], work:Path, tag:str):
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
    raise RuntimeError(f'cannot resolve required outputs from {x}')

def parse_fasta(p):
    d={}; h=None; seq=[]
    with open(p) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: d[h]=''.join(seq)
                h=line[1:].split()[0]; seq=[]
            else: seq.append(line)
    if h is not None: d[h]=''.join(seq)
    return d

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('f1_input'); ap.add_argument('f1r1_input'); ap.add_argument('--outdir',default='outputs'); ap.add_argument('--config',default='config/FROZEN_REFERENCE_LOCKS.tsv')
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); work=Path('work'); work.mkdir(exist_ok=True)
    f1_names=['STAGE05F1_RESCUE_PRIORITY.tsv','STAGE05F1_SIXFRAME_PEPTIDE_INSTANCES.tsv','STAGE05F1_UNIQUE_SIXFRAME_PEPTIDES.faa','STAGE05F1_GATE.json']
    r1_names=['STAGE05F1R1_STRICT_FORMAL_GENE_MODEL_RESCUES.tsv','STAGE05F1R1_PRIORITY_INSTANCE_AUDIT.tsv','STAGE05F1R1_GATE.json']
    f1,route1=resolve_dir_or_tar(Path(a.f1_input).resolve(),f1_names,work,'f1')
    r1,route2=resolve_dir_or_tar(Path(a.f1r1_input).resolve(),r1_names,work,'r1')
    g1=json.loads((f1/'STAGE05F1_GATE.json').read_text()); g2=json.loads((r1/'STAGE05F1R1_GATE.json').read_text())
    if not g1.get('gate_pass') or not g2.get('gate_pass'): raise RuntimeError('upstream gate not PASS')
    if int(g2.get('strict_formal_gene_model_rescue_unique_aa_n',-1))!=39: raise RuntimeError('expected frozen strict formal count 39')
    locks=read_tsv(Path(a.config)); cap=set(); rtacc=set(); d1acc=set(); label={}
    for r in locks:
        t=r['lock_type']; k=r['key']; label[(t,k)]=r.get('label','')
        if t=='STAGE05C1_CAPTURED_RT_CANDIDATE': cap.add(k)
        elif t=='STAGE05C1_RT_CORE_ACCESSION': rtacc.add(k)
        elif t=='PRIOR_D1_FORMAL_ACCESSION': d1acc.add(k)
    strict=read_tsv(r1/'STAGE05F1R1_STRICT_FORMAL_GENE_MODEL_RESCUES.tsv')
    inst=read_tsv(r1/'STAGE05F1R1_PRIORITY_INSTANCE_AUDIT.tsv')
    seqs=parse_fasta(f1/'STAGE05F1_UNIQUE_SIXFRAME_PEPTIDES.faa')
    # query_id->sequence; exact hash sanity
    for q,s in seqs.items():
        pass
    rows=[]
    for r in strict:
        cands=splitsemi(r['candidate_keys']); accs=splitsemi(r['accessions'])
        exact_caps=sorted(set(cands)&cap); rt_same=sorted(set(accs)&rtacc); d1_same=sorted(set(accs)&d1acc)
        if exact_caps:
            provenance='KNOWN_CAPTURED_RT_CALIBRATION'
        elif d1_same:
            provenance='PRIOR_D1_FORMAL_CONTEXT'
        elif rt_same:
            provenance='SAME_GENOME_AS_FROZEN_RT_CORE_OTHER_LOCUS'
        else:
            provenance='NOVELTY_FIRST_STRICT_RESCUE'
        highcov=(flt(r.get('query_coverage'))>=0.50 and flt(r.get('model_coverage'))>=0.25)
        yxdd=asbool(r.get('yxdd_like')); replicated=int(r.get('genome_n') or 0)>=2
        pure_no=(r.get('strict_aa_status')=='PURE_NO_SAME_STRAND_OVERLAP_NONEDGE')
        if provenance=='KNOWN_CAPTURED_RT_CALIBRATION': queue='CALIBRATION_ONLY'
        elif provenance=='PRIOR_D1_FORMAL_CONTEXT': queue='PRIOR_REVIEW_BRIDGE'
        elif yxdd and pure_no: queue='P1_YXDD_PURE_MISSED_ORF'
        elif replicated: queue='P2_CROSS_GENOME_REPLICATED'
        elif highcov: queue='P3_HIGHCOV_SINGLETON'
        else: queue='P4_OTHER_STRICT_FORMAL'
        q=r['query_id']; s=seqs.get(q,'')
        if not s: raise RuntimeError(f'missing sequence for {q}')
        if hashlib.sha256(s.encode()).hexdigest()!=r['aa_sha256']: raise RuntimeError(f'hash mismatch {q}')
        x=dict(r); x.update({'f2a_provenance_class':provenance,'exact_known_rt_candidate_keys':';'.join(exact_caps),'same_genome_frozen_rt_core_accessions':';'.join(rt_same),'prior_d1_formal_accessions':';'.join(d1_same),'highcov_descriptive':str(highcov),'replicated_multi_genome':str(replicated),'pure_no_same_strand_overlap':str(pure_no),'compute_queue':queue})
        rows.append(x)
    order={'CALIBRATION_ONLY':0,'PRIOR_REVIEW_BRIDGE':1,'P1_YXDD_PURE_MISSED_ORF':2,'P2_CROSS_GENOME_REPLICATED':3,'P3_HIGHCOV_SINGLETON':4,'P4_OTHER_STRICT_FORMAL':5}
    rows.sort(key=lambda r:(order[r['compute_queue']],flt(r.get('best_full_evalue'),9e99),-flt(r.get('model_coverage')),r['query_id']))
    write_tsv(out/'STAGE05F2A_STRICT_FORMAL_CALIBRATED.tsv',rows)
    known=[r for r in rows if r['f2a_provenance_class']=='KNOWN_CAPTURED_RT_CALIBRATION']
    novelty=[r for r in rows if r['f2a_provenance_class']=='NOVELTY_FIRST_STRICT_RESCUE']
    prior=[r for r in rows if r['f2a_provenance_class']!='NOVELTY_FIRST_STRICT_RESCUE' and r['f2a_provenance_class']!='KNOWN_CAPTURED_RT_CALIBRATION']
    write_tsv(out/'STAGE05F2A_KNOWN_RT_CALIBRATION.tsv',known)
    write_tsv(out/'STAGE05F2A_NOVELTY_FIRST_STRICT_RESCUES.tsv',novelty)
    write_tsv(out/'STAGE05F2A_PRIOR_CONTEXT_BRIDGES.tsv',prior)
    compute=[r for r in rows if r['compute_queue'] in ('P1_YXDD_PURE_MISSED_ORF','P2_CROSS_GENOME_REPLICATED','P3_HIGHCOV_SINGLETON')]
    write_tsv(out/'STAGE05F2A_COMPUTE_PRIORITY.tsv',compute)
    # instance-level records for novelty rows
    novelty_hash={r['aa_sha256'] for r in novelty}; inst_out=[r for r in inst if r.get('aa_sha256') in novelty_hash]
    write_tsv(out/'STAGE05F2A_NOVELTY_FIRST_INSTANCE_MAP.tsv',inst_out)
    # fasta outputs
    for fn,subset in [('STAGE05F2A_ALL_STRICT_FORMAL.faa',rows),('STAGE05F2A_NOVELTY_FIRST_STRICT.faa',novelty),('STAGE05F2A_COMPUTE_PRIORITY.faa',compute)]:
        with (out/fn).open('w') as f:
            for r in subset:
                s=seqs[r['query_id']]
                f.write(f">{r['query_id']} sha256={r['aa_sha256']} len={len(s)} queue={r['compute_queue']} profile={r['best_rvt_target']} provenance={r['f2a_provenance_class']}\n{s}\n")
    # calibration distribution
    calib={'known_control_n':len(known),'known_control_yxdd_n':sum(asbool(r['yxdd_like']) for r in known),'known_control_full_evalue_max':max([flt(r['best_full_evalue']) for r in known],default=None),'known_control_model_coverage_min':min([flt(r['model_coverage']) for r in known],default=None),'known_control_query_coverage_min':min([flt(r['query_coverage']) for r in known],default=None)}
    (out/'STAGE05F2A_KNOWN_CONTROL_CALIBRATION.json').write_text(json.dumps(calib,indent=2)+'\n')
    counts=Counter(r['f2a_provenance_class'] for r in rows); queues=Counter(r['compute_queue'] for r in rows)
    summary={'stage':'Stage05F2A','version':'v1.0.0','gate':'PASS_STAGE05F2A_STRICT_RESCUE_CALIBRATION_AND_COMPUTE_SHORTLIST','gate_pass':True,'f1_route':route1,'f1r1_route':route2,'strict_formal_unique_aa_n':len(rows),'provenance_class_counts':dict(counts),'compute_queue_counts':dict(queues),'known_calibration_n':len(known),'novelty_first_strict_n':len(novelty),'compute_priority_n':len(compute),'reference_lock_candidate_n':len(cap),'reference_lock_rt_accession_n':len(rtacc),'prior_d1_accession_n':len(d1acc),'new_hmm_search_performed':False,'new_structure_prediction_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'candidate_membership_changed':False,'claim_boundary':'Post-rescue calibration and compute triage only. Novelty-first strict rescues remain RT-like candidates, not biochemical RTs or ART systems. Known/near-known RT contexts are retained as positive calibration, not discoveries.'}
    (out/'STAGE05F2A_GATE.json').write_text(json.dumps(summary,indent=2)+'\n')
    notes=['# Stage05F2A calibrated rescue interpretation\n\n',f'- Strict formal non-edge rescue unique AA: {len(rows)}.\n',f'- Known captured RT calibration: {len(known)}.\n',f'- Novelty-first strict rescue: {len(novelty)}.\n',f'- Descriptive compute-priority subset: {len(compute)}.\n\n','Known ART/RT-related loci recovered by six-frame rescue are calibration positives and must not be counted as novel rescue discoveries. Prior Stage05D1 formal contexts are separately bridged. The remaining novelty-first strict set is retained without deleting weak members; compute queues only prioritize the next orthogonal analyses.\n']
    (out/'STAGE05F2A_MANUSCRIPT_NOTES.md').write_text(''.join(notes))
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
