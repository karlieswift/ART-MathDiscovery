#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,gzip,json,tarfile
from collections import defaultdict,Counter
from pathlib import Path
FORMAL='FORMAL_PROFILE_GENE_CALLER_RESCUE'
EXPL='EXPLORATORY_PROFILE_GENE_CALLER_RESCUE'
REP='HMM_UNDER_DETECTED_MOTIF_REPLICATED_RESCUE_RESERVOIR'
SING='HMM_UNDER_DETECTED_MOTIF_SINGLETON_RESCUE_RESERVOIR'

def open_text(p): return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p,'r',newline='')
def read_tsv(p):
    with open_text(p) as f: return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,cols=None):
    if cols is None: cols=list(rows[0]) if rows else []
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
def resolve_f1(x:Path, work:Path):
    if x.is_dir():
        out=x/'outputs'
        if (out/'STAGE05F1_RESCUE_PRIORITY.tsv').exists(): return out,'live_directory'
    if x.is_file() and str(x).endswith('.tar.gz'):
        dest=work/'f1_extract'; dest.mkdir(parents=True,exist_ok=True)
        names={'STAGE05F1_RESCUE_PRIORITY.tsv','STAGE05F1_SIXFRAME_PEPTIDE_INSTANCES.tsv','STAGE05F1_RVT_ALL.domtblout','STAGE05F1_GATE.json'}
        with tarfile.open(x,'r:gz') as tf:
            for m in tf.getmembers():
                if m.isfile() and Path(m.name).name in names:
                    fh=tf.extractfile(m)
                    if fh: (dest/Path(m.name).name).write_bytes(fh.read())
        if (dest/'STAGE05F1_RESCUE_PRIORITY.tsv').exists(): return dest,'result_tar'
    raise RuntimeError(f'cannot resolve Stage05F1 outputs from {x}')
def parse_domtbl(p):
    byq=defaultdict(list)
    if not p.exists(): return byq
    with open(p) as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            x=line.rstrip('\n').split(maxsplit=22)
            if len(x)<22: continue
            try:
                r={'target_name':x[0],'target_hmm_len':int(x[2]),'query_id':x[3],'query_len_domtbl':int(x[5]),'full_evalue':float(x[6]),'full_score':float(x[7]),'domain_c_evalue':float(x[11]),'domain_i_evalue':float(x[12]),'domain_score':float(x[13]),'hmm_from':int(x[15]),'hmm_to':int(x[16]),'ali_from':int(x[17]),'ali_to':int(x[18]),'env_from':int(x[19]),'env_to':int(x[20]),'domain_acc':float(x[21]),'description':x[22] if len(x)>22 else ''}
                r['query_coverage']=(r['ali_to']-r['ali_from']+1)/max(1,r['query_len_domtbl'])
                r['model_coverage']=(r['hmm_to']-r['hmm_from']+1)/max(1,r['target_hmm_len'])
                byq[r['query_id']].append(r)
            except Exception: pass
    for q in byq: byq[q].sort(key=lambda r:(r['full_evalue'],r['domain_i_evalue'],-r['domain_score']))
    return byq
def splitsemi(x): return [z for z in str(x or '').split(';') if z]
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('f1_input'); ap.add_argument('project_root',nargs='?',default='.'); ap.add_argument('--outdir',default='outputs')
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); work=Path('work'); work.mkdir(exist_ok=True)
    src,route=resolve_f1(Path(a.f1_input).resolve(),work)
    pri=read_tsv(src/'STAGE05F1_RESCUE_PRIORITY.tsv'); inst=read_tsv(src/'STAGE05F1_SIXFRAME_PEPTIDE_INSTANCES.tsv'); dom=parse_domtbl(src/'STAGE05F1_RVT_ALL.domtblout')
    gate=json.loads((src/'STAGE05F1_GATE.json').read_text()) if (src/'STAGE05F1_GATE.json').exists() else {}
    pbyhash={r['aa_sha256']:r for r in pri}; inst_by_hash=defaultdict(list)
    for r in inst:
        if r.get('aa_sha256') in pbyhash: inst_by_hash[r['aa_sha256']].append(r)
    aa_summary=[]; audited=[]; class_order=[FORMAL,EXPL,REP,SING]
    for h,p in pbyhash.items():
        cls=p.get('rescue_class','')
        if cls not in class_order: continue
        q=p.get('query_id',''); best=dom.get(q,[{}])[0] if dom.get(q) else {}; rr=inst_by_hash.get(h,[])
        nonedge_no=nonedge_boundary=edge_no=edge_boundary=contained=0; accs=set(); cands=set()
        for x in rr:
            rel=x.get('prodigal_relation',''); edge=x.get('edge_partial','') in ('True','true','1','TRUE')
            if rel=='NO_SAME_STRAND_PRODIGAL_OVERLAP':
                if edge: edge_no+=1
                else: nonedge_no+=1
            elif rel=='BOUNDARY_OR_PARTIAL_DISCORDANT_WITH_PRODIGAL':
                if edge: edge_boundary+=1
                else: nonedge_boundary+=1
            elif rel=='CONTAINED_IN_EXISTING_PRODIGAL_ORF': contained+=1
            accs.add(x.get('sequence_id','')); cands.update(splitsemi(x.get('candidate_keys','')))
            audited.append({'rescue_class':cls,'query_id':q,'aa_sha256':h,'aa_length':p.get('aa_length',''),'sequence_id':x.get('sequence_id',''),'start0':x.get('start0',''),'end0':x.get('end0',''),'strand':x.get('strand',''),'frame':x.get('frame',''),'edge_partial':x.get('edge_partial',''),'prodigal_relation':rel,'candidate_keys':x.get('candidate_keys',''),'yxdd_like':p.get('yxdd_like',''),'best_rvt_target':best.get('target_name',p.get('best_rvt_target','')),'best_full_evalue':best.get('full_evalue',p.get('best_full_evalue','')),'best_domain_i_evalue':best.get('domain_i_evalue',p.get('best_domain_i_evalue','')),'query_coverage':best.get('query_coverage',''),'model_coverage':best.get('model_coverage',''),'domain_score':best.get('domain_score',''),'domain_acc':best.get('domain_acc','')})
        rescue_total=nonedge_no+nonedge_boundary+edge_no+edge_boundary
        if rr and nonedge_no==len(rr): strict='PURE_NO_SAME_STRAND_OVERLAP_NONEDGE'
        elif rr and contained==0 and rescue_total==len(rr) and (nonedge_no+nonedge_boundary)>0: strict='PURE_GENE_MODEL_RESCUE_WITH_NONEDGE'
        elif contained>0 and rescue_total>0 and (nonedge_no+nonedge_boundary)>0: strict='MIXED_FANOUT_EXISTING_AND_RESCUED_NONEDGE'
        elif rr and (nonedge_no+nonedge_boundary)==0 and (edge_no+edge_boundary)>0: strict='EDGE_PARTIAL_ONLY_RESCUE'
        elif rr and contained==len(rr): strict='EXISTING_ORF_ONLY'
        else: strict='OTHER_MIXED'
        aa_summary.append({'rescue_class':cls,'strict_aa_status':strict,'query_id':q,'aa_sha256':h,'aa_length':p.get('aa_length',''),'instance_n':len(rr),'genome_n':len([x for x in accs if x]),'candidate_n':len(cands),'nonedge_no_same_strand_overlap_n':nonedge_no,'nonedge_boundary_discordant_n':nonedge_boundary,'edge_no_same_strand_overlap_n':edge_no,'edge_boundary_discordant_n':edge_boundary,'contained_existing_orf_n':contained,'yxdd_like':p.get('yxdd_like',''),'best_rvt_target':best.get('target_name',p.get('best_rvt_target','')),'best_full_evalue':best.get('full_evalue',p.get('best_full_evalue','')),'best_domain_i_evalue':best.get('domain_i_evalue',p.get('best_domain_i_evalue','')),'query_coverage':best.get('query_coverage',''),'model_coverage':best.get('model_coverage',''),'domain_score':best.get('domain_score',''),'domain_acc':best.get('domain_acc',''),'candidate_keys':p.get('candidate_keys',''),'accessions':p.get('accessions','')})
    def fev(r):
        try: return float(r['best_full_evalue'])
        except: return 9e99
    aa_summary.sort(key=lambda r:(class_order.index(r['rescue_class']),fev(r),-float(r['query_coverage'] or 0),-float(r['model_coverage'] or 0)))
    write_tsv(out/'STAGE05F1R1_PRIORITY_AA_STRICT_AUDIT.tsv',aa_summary); write_tsv(out/'STAGE05F1R1_PRIORITY_INSTANCE_AUDIT.tsv',audited)
    formal=[r for r in aa_summary if r['rescue_class']==FORMAL]
    strict=[r for r in formal if r['strict_aa_status'] in ('PURE_NO_SAME_STRAND_OVERLAP_NONEDGE','PURE_GENE_MODEL_RESCUE_WITH_NONEDGE')]
    mixed=[r for r in formal if r['strict_aa_status']=='MIXED_FANOUT_EXISTING_AND_RESCUED_NONEDGE']; edge=[r for r in formal if r['strict_aa_status']=='EDGE_PARTIAL_ONLY_RESCUE']
    highcov=[r for r in strict if float(r['query_coverage'] or 0)>=0.5 and float(r['model_coverage'] or 0)>=0.25]
    expl=[r for r in aa_summary if r['rescue_class']==EXPL and (int(r['nonedge_no_same_strand_overlap_n'])+int(r['nonedge_boundary_discordant_n'])>0)]
    rep=[r for r in aa_summary if r['rescue_class']==REP and (int(r['nonedge_no_same_strand_overlap_n'])+int(r['nonedge_boundary_discordant_n'])>0)]
    write_tsv(out/'STAGE05F1R1_STRICT_FORMAL_GENE_MODEL_RESCUES.tsv',strict); write_tsv(out/'STAGE05F1R1_STRICT_FORMAL_HIGHCOV_DESCRIPTIVE.tsv',highcov); write_tsv(out/'STAGE05F1R1_EXPLORATORY_PROFILE_NONEDGE.tsv',expl); write_tsv(out/'STAGE05F1R1_HMM_UNDERDETECTED_REPLICATED_NONEDGE.tsv',rep)
    prof=Counter(r['best_rvt_target'] for r in formal if r['best_rvt_target']); write_tsv(out/'STAGE05F1R1_FORMAL_PROFILE_DISTRIBUTION.tsv',[{'rvt_profile':k,'formal_unique_aa_n':v} for k,v in sorted(prof.items(),key=lambda kv:(-kv[1],kv[0]))],['rvt_profile','formal_unique_aa_n'])
    strict_g=set(); strict_c=set()
    for r in strict: strict_g.update(splitsemi(r['accessions'])); strict_c.update(splitsemi(r['candidate_keys']))
    summary={'stage':'Stage05F1R1','version':'v1.0.0','gate':'PASS_STAGE05F1R1_STRICT_INSTANCE_DECONVOLUTION','gate_pass':True,'f1_input_route':route,'f1_gate':gate.get('gate',''),'f1_gate_pass':gate.get('gate_pass',None),'f1_formal_gene_caller_rescue_unique_aa_n':len(formal),'formal_unique_aa_strict_status_counts':dict(Counter(r['strict_aa_status'] for r in formal)),'strict_formal_gene_model_rescue_unique_aa_n':len(strict),'strict_formal_gene_model_rescue_genome_n':len(strict_g),'strict_formal_gene_model_rescue_candidate_n':len(strict_c),'mixed_fanout_formal_unique_aa_n':len(mixed),'edge_partial_only_formal_unique_aa_n':len(edge),'strict_formal_highcov_descriptive_unique_aa_n':len(highcov),'highcov_descriptive_rule':'query_coverage>=0.50 and model_coverage>=0.25; descriptive only, not a frozen biological threshold','exploratory_profile_nonedge_unique_aa_n':len(expl),'hmm_underdetected_replicated_nonedge_unique_aa_n':len(rep),'new_hmm_search_performed':False,'new_sixframe_translation_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'candidate_membership_changed':False,'claim_boundary':'Strict deconvolution of Stage05F1 evidence. Even strict formal rescue remains a candidate missed/boundary-missed RT-like peptide, not biochemical RT or ART proof.'}
    (out/'STAGE05F1R1_GATE.json').write_text(json.dumps(summary,indent=2)+'\n')
    note=(f'# Stage05F1R1 strict instance audit\n\nStage05F1 reported {len(formal)} unique-AA formal profile gene-caller rescue classes. F1R1 decomposes those classes across all peptide instances.\n\n'
          f'- strict formal gene-model rescue unique AA: {len(strict)}\n'
          f'- strict formal rescue genomes: {len(strict_g)}\n'
          f'- strict formal rescue candidate neighborhoods: {len(strict_c)}\n'
          f'- mixed existing/rescued fanout AA: {len(mixed)}\n'
          f'- edge-partial-only formal AA: {len(edge)}\n'
          f'- descriptive high-coverage strict formal AA: {len(highcov)}\n'
          f'- exploratory-profile nonedge AA: {len(expl)}\n'
          f'- HMM-under-detected replicated nonedge AA: {len(rep)}\n\n'
          'Do not compare the F1 formal count directly with the five Stage05D1 hits: the universes and translation rules differ.\n')
    (out/'STAGE05F1R1_MANUSCRIPT_NOTES.md').write_text(note)
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
