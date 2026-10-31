#!/usr/bin/env python3
import csv, json, sys, hashlib
from pathlib import Path
from common import *

def main():
    if len(sys.argv)!=6:
        raise SystemExit('usage: prepare.py B4D_DIR OUTDIR INPUTDIR FROZENDIR LOCK')
    b4d,outdir,inputdir,frozendir,lockp=map(Path,sys.argv[1:])
    outdir.mkdir(parents=True,exist_ok=True); inputdir.mkdir(parents=True,exist_ok=True); frozendir.mkdir(parents=True,exist_ok=True)
    board=b4d/'outputs/F3B4D_FINAL_SEVEN_FAMILY_COMPUTATIONAL_EVIDENCE_BOARD.tsv'
    reps=b4d/'outputs/F3B4D_ARCHITECTURE_REPRESENTATIVE_PER_FAMILY_7.tsv'
    rtfa=b4d/'inputs/F3B4D_RT_ARCHITECTURE_REPRESENTATIVES_7.faa'
    pfa=b4d/'inputs/F3B4D_PARTNER_REPRESENTATIVES_7.faa'
    cfa=b4d/'inputs/F3B4D_CORE_ARCHITECTURE_CONTEXT_16.fna'
    gate=b4d/'outputs/STAGE05F3B4D_GATE.json'
    for p in [board,reps,rtfa,pfa,cfa,gate]:
        if not p.exists(): raise SystemExit(f'[FATAL] missing {p}')
    g=json.load(open(gate)); assert g['gate_pass'] is True
    br=read_tsv(board)
    targets=[r for r in br if r['structure_tier'] in ('STRONG','SUGGESTIVE') and r['computational_art_like_status'].startswith('RT_FOLD_SUPPORTED_STRICT_ARRAY_RECURRENT_PARTNER')]
    targets=sorted(targets,key=lambda r:int(r['computational_priority_rank']))
    if len(targets)!=5: raise SystemExit(f'[FATAL] expected 5 RT-fold-supported ART-like targets, got {len(targets)}')
    fams=[r['family_rep_id'] for r in targets]; famset=set(fams)
    rr={r['family_rep_id']:r for r in read_tsv(reps) if r['family_rep_id'] in famset}
    if set(rr)!=famset: raise SystemExit('[FATAL] representative table missing target family')
    rt,rt_desc=read_fasta(rtfa); pp,p_desc=read_fasta(pfa); cc,c_desc=read_fasta(cfa)
    # map family from descriptions for protein reps
    rt_byfam={family_from_desc(rt_desc[k]):(k,rt[k],rt_desc[k]) for k in rt if family_from_desc(rt_desc[k]) in famset}
    p_byfam={family_from_desc(p_desc[k]):(k,pp[k],p_desc[k]) for k in pp if family_from_desc(p_desc[k]) in famset}
    if set(rt_byfam)!=famset or set(p_byfam)!=famset: raise SystemExit('[FATAL] RT/partner FASTA family mapping incomplete')
    # core architecture representative candidate IDs are frozen in reps table.
    arch=[]; rt5=[]; p5=[]; rows=[]
    for rank,fam in enumerate(fams,1):
        rep=rr[fam]; cid=rep['candidate_id']
        if cid not in cc: raise SystemExit(f'[FATAL] representative core DNA missing {cid}')
        rid,rseq,_=rt_byfam[fam]; pid,pseq,_=p_byfam[fam]
        if rid!=cid: raise SystemExit(f'[FATAL] RT representative mismatch family={fam}: fasta={rid} table={cid}')
        rt5.append(f'>{cid} family={fam} role=RT_FOLD_SUPPORTED_ART_LIKE\n{wrap(rseq)}')
        p5.append(f'>{pid} family={fam} role=RECURRENT_DOWNSTREAM_PARTNER_CANDIDATE\n{wrap(pseq)}')
        arch.append(f'>{cid} family={fam} role=STRICT_ARRAY_RT_PARTNER_ARCHITECTURE_CORE\n{wrap(cc[cid])}')
        rows.append({
          'priority_rank':rank,'family_rep_id':fam,'rt_candidate_id':cid,'rt_aa_length':len(rseq),'rt_sha256':hashlib.sha256(rseq.encode()).hexdigest(),
          'partner_orf_id':pid,'partner_aa_length':len(pseq),'partner_sha256':hashlib.sha256(pseq.encode()).hexdigest(),
          'architecture_core_bp':len(cc[cid]),'architecture_core_sha256':hashlib.sha256(cc[cid].encode()).hexdigest(),
          'structure_tier':targets[rank-1]['structure_tier'],'source_accession':rep['accession'],
          'repeat_motif':targets[rank-1]['modal_exact_motif_core'],'spacer_bp_median':targets[rank-1]['spacer_bp_median'],
          'partner_phrog':targets[rank-1]['partner_phrog_top_target'],'partner_pfam':targets[rank-1]['partner_pfam_top_target']
        })
    (inputdir/'F3B4E_RT5.faa').write_text(''.join(rt5),encoding='ascii')
    (inputdir/'F3B4E_PARTNER5.faa').write_text(''.join(p5),encoding='ascii')
    (inputdir/'F3B4E_ARCHITECTURE5.fna').write_text(''.join(arch),encoding='ascii')
    write_tsv(outdir/'F3B4E_TARGETS_5.tsv',rows)
    art8=[
      ('ART_21','SA1','MW218148.1','author Methods array case; Figure 3 ART_21'),
      ('','MarsHill','MW248466.1','author Methods array case; initial jumbo-phage context'),
      ('ART_19','Madawaska','MW349129.1','Figure 3 ART_19'),
      ('ART_20','LY01','OR836606.1','Figure 3 ART_20'),
      ('','S6','LC680885.1','author Methods array case'),
      ('','PALS_2','MN091626.1','author Methods array case'),
      ('','UFV_DC4','MZ779063.1','author Methods array case'),
      ('ART_17','PB50','OR770614.1','Figure 3 ART_17')]
    arows=[{'art_id':a,'source_name':n,'accession':x,'public_evidence':e} for a,n,x,e in art8]
    write_tsv(inputdir/'PUBLIC_ART8_MANIFEST.tsv',arows)
    (inputdir/'PUBLIC_ART8_ACCESSIONS.txt').write_text('\n'.join(x for _,_,x,_ in art8)+'\n')
    lock={
      'stage':'Stage05F3B4E-implementation-freeze','gate':'PASS_F3B4E_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_RULES_FROZEN_PRE_RESULTS','gate_pass':True,
      'target_rule':'all and only five F3B4D final families with frozen STRONG/SUGGESTIVE RT fold plus strict upstream array plus recurrent downstream partner candidate','target_family_n':5,'target_family_ids':fams,
      'gap1_global_public_sequence_audit':'NCBI remote BLASTP against nr for RT and partner representatives plus remote MEGABLAST against nt for frozen architecture cores; query/result date is recorded at runtime',
      'protein_exact_rule':'pident>=99.999%, qcov>=99%, tcov>=99%','protein_high_identity_near_full_rule':'pident>=90%, qcov>=80%, tcov>=80%; same rule as F3B3B local known-system audit',
      'architecture_exact_like_rule':'pident>=99.9% and qcov>=95% in one HSP','architecture_near_full_rule':'pident>=95% and qcov>=80% in one HSP',
      'gap2_public_art_crosswalk':'exact public named ART8 genome panel; local TBLASTN of candidate RT and partner proteins and DC-MEGABLAST/BLASTN of architecture core',
      'anthropic_art_cluster_threshold_context':'Anthropic reports 95 ART RT clusters at 90% identity. Because the public ART8 crosswalk searches whole nucleotide genomes by TBLASTN rather than a complete ART protein catalogue, we report 90%-identity/80%-query-coverage compatibility only; it is not a formal ART-cluster assignment.',
      'public_art8_accession_n':8,'full_art95_machine_readable_catalogue_present':False,
      'full_art95_crosswalk_claim_authorized':False,'anthropic_unreported_art_claim_authorized':False,
      'absolute_global_novelty_claim_authorized_from_ncbi_nr_nt_alone':False,
      'b4d_gate_sha256':sha256_file(gate),'b4d_final_board_sha256':sha256_file(board),'b4d_rep_table_sha256':sha256_file(reps),
      'query_rt5_sha256':sha256_file(inputdir/'F3B4E_RT5.faa'),'query_partner5_sha256':sha256_file(inputdir/'F3B4E_PARTNER5.faa'),'query_architecture5_sha256':sha256_file(inputdir/'F3B4E_ARCHITECTURE5.fna'),
      'public_art8_manifest_sha256':sha256_file(inputdir/'PUBLIC_ART8_MANIFEST.tsv'),
      'claim_boundary':'This closes a dated public-NCBI sequence/locus audit and the exact public named ART8 subset crosswalk. It cannot by itself prove absolute global novelty or absence from the unpublished/non-machine-readable remainder of ART95.'}
    Path(lockp).write_text(json.dumps(lock,indent=2)+'\n')
    print(json.dumps(lock,indent=2))
    print('\n=== TARGETS 5 ==='); print((outdir/'F3B4E_TARGETS_5.tsv').read_text().rstrip())
    print('PASS_F3B4E_PREPARE_COMPLETE')

if __name__=='__main__': main()
