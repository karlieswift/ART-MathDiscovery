#!/usr/bin/env python3
import csv, json, hashlib, importlib.util, sys
from pathlib import Path
from collections import defaultdict


def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def loadmod(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def read_fasta(path):
    out={}; h=None; seq=[]
    with open(path,encoding='utf-8',errors='ignore') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: out[h]=''.join(seq).upper()
                h=line[1:].split()[0]; seq=[]
            else: seq.append(line)
        if h is not None: out[h]=''.join(seq).upper()
    return out

def wrap(s,w=80):
    return '\n'.join(s[i:i+w] for i in range(0,len(s),w))+'\n'

def boolv(x): return str(x).lower()=='true'

def main():
    if len(sys.argv)!=15:
        raise SystemExit('usage: prepare_targets.py B4C_OUT R4_INSTANCE_MAP ORFDB ORFDB_INDEX GENOMEDB GENOMEDB_INDEX ER1_SUMMARY VENDOR_HELPER OUTDIR INPUTDIR FROZENDIR STRUCT_YAML_DIR LOCK_JSON ROOT')
    b4c_out, imap_p, orfdb, orfidx, gdb, gidx, er1_sum, vendor, outdir, inputdir, frozendir, yaml_dir, lock_json, root = sys.argv[1:]
    b4c=Path(b4c_out); out=Path(outdir); inp=Path(inputdir); frz=Path(frozendir); ydir=Path(yaml_dir)
    for d in (out,inp,frz,ydir): d.mkdir(parents=True,exist_ok=True)
    v=loadmod('ctx',vendor)

    famsum_p=b4c/'F3B4C2_RT_FAMILY_REPEAT_PARTNER_ARCHITECTURE_SUMMARY_4912.tsv'
    arch_p=b4c/'F3B4C2_RT_REPEAT_PARTNER_LOCUS_ARCHITECTURE_10598.tsv'
    arr_p=b4c/'F3B4C1_ALL_DETECTED_REPEAT_ARRAYS.tsv'
    loc_p=b4c/'F3B4C1_GLOBAL_FORMAL_RT_LOCUS_REPEAT_CENSUS_10598.tsv'
    near_p=b4c/'F3B4C2_NEAREST_DOWNSTREAM_ORF_PER_RT_LOCUS.tsv'
    pfa_p=b4c/'F3B4C2_NEAREST_DOWNSTREAM_ORF_PROTEINS.faa'
    gate_p=b4c/'STAGE05F3B4C2_GATE.json'
    for p in [famsum_p,arch_p,arr_p,loc_p,near_p,pfa_p,gate_p,Path(imap_p),Path(orfdb),Path(orfidx),Path(gdb),Path(gidx),Path(er1_sum),Path(vendor)]:
        if not p.exists(): raise SystemExit(f'FAIL missing {p}')
    gate=json.load(open(gate_p)); assert gate['gate_pass'] is True
    assert gate['family_n_with_strict_repeat_and_recurrent_nearest_partner_cluster']==7

    famsum=list(csv.DictReader(open(famsum_p,encoding='utf-8'),delimiter='\t'))
    targets=[r for r in famsum if boolv(r['strict_array_plus_recurrent_partner_cluster'])]
    if len(targets)!=7: raise SystemExit(f'FAIL target family n={len(targets)}')
    targets=sorted(targets,key=lambda r:(-int(r['strict_array_dominant_partner_cluster_accession_n']),-int(r['strict_array_locus_n']),r['family_rep_id']))
    famset={r['family_rep_id'] for r in targets}
    dominant_cluster={r['family_rep_id']:r['strict_array_dominant_partner_cluster_id'] for r in targets}

    loc={r['candidate_id']:r for r in csv.DictReader(open(loc_p,encoding='utf-8'),delimiter='\t') if r['family_rep_id'] in famset}
    arch={r['candidate_id']:r for r in csv.DictReader(open(arch_p,encoding='utf-8'),delimiter='\t') if r['family_rep_id'] in famset}
    near={r['candidate_id']:r for r in csv.DictReader(open(near_p,encoding='utf-8'),delimiter='\t') if r['family_rep_id'] in famset}
    arrays=[r for r in csv.DictReader(open(arr_p,encoding='utf-8'),delimiter='\t') if r['family_rep_id'] in famset and boolv(r['strict_published_art_geometry'])]
    arr_by=defaultdict(list)
    for r in arrays: arr_by[r['candidate_id']].append(r)

    core=[]
    for cid,a in arch.items():
        fam=a['family_rep_id']
        if boolv(a['has_strict_published_art_geometry']) and a['nearest_downstream_partner_cluster_id']==dominant_cluster[fam] and a['nearest_downstream_orf_id']:
            core.append(cid)
    core=sorted(core)
    if len(core)!=16: raise SystemExit(f'FAIL core recurrent strict locus n={len(core)} expected 16')
    if any(cid not in arr_by for cid in core): raise SystemExit('FAIL core locus lacks strict array record')

    # instance provenance for RT peptide extraction
    want=set(core)
    # also need all family representative IDs in case present as candidate instances, but architecture reps are chosen from core loci below.
    im={}
    with open(imap_p,encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['instance_id'] in want: im[r['instance_id']]=r
    if set(core)!=set(im): raise SystemExit(f'FAIL instance map missing {sorted(set(core)-set(im))}')
    oi=v.load_index(orfidx,{im[c]['orfdb_key'] for c in core})
    gi=v.load_index(gidx,{im[c]['parent_contig_key'] for c in core})
    peps={}; contigs={}
    with open(orfdb,'rb',buffering=0) as f:
        for k,(o,l) in oi.items(): peps[k]=v.fetch_db(f,o,l).rstrip('*')
    with open(gdb,'rb',buffering=0) as f:
        for k,(o,l) in gi.items(): contigs[k]=v.fetch_db(f,o,l)
    partner_fa=read_fasta(pfa_p)

    context_rows=[]; core_fa=[]; padded_fa=[]; rt_core_fa=[]; partner_core_fa=[]
    for cid in core:
        l=loc[cid]; a=arch[cid]; n=near[cid]; fam=l['family_rep_id']; p=im[cid]
        seq=contigs[p['parent_contig_key']]; pep=peps[p['orfdb_key']]
        s=int(l['candidate_start0']); e=int(l['candidate_end0']); strand=l['strand']
        if (e-s)!=3*len(pep): raise SystemExit(f'FAIL RT nt/aa length mismatch {cid}')
        if strand=='+':
            up=seq[max(0,s-10000):s]; rt=seq[s:e]; dn=seq[e:min(len(seq),e+10000)]
        elif strand=='-':
            up=v.revcomp(seq[e:min(len(seq),e+10000)]); rt=v.revcomp(seq[s:e]); dn=v.revcomp(seq[max(0,s-10000):s])
        else: raise SystemExit(f'FAIL strand {cid}')
        comb=up+rt+dn; off=len(up); rtnt=len(rt)
        relstarts=[]
        for ar in arr_by[cid]:
            relstarts += [int(x) for x in ar['candidate_relative_repeat_starts'].split(',') if x!='']
        left_rel=min(relstarts)
        base=3 if boolv(n['rt_terminal_stop_present']) else 0
        if not n['local_end0_including_stop_after_rt_stop']: raise SystemExit(f'FAIL partner end missing {cid}')
        partner_end_rel=rtnt+base+int(n['local_end0_including_stop_after_rt_stop'])
        core_s=off+left_rel; core_e=off+partner_end_rel
        if core_s<0 or core_e>len(comb) or core_s>=core_e: raise SystemExit(f'FAIL core bounds {cid} {core_s} {core_e} len={len(comb)}')
        core_seq=comb[core_s:core_e]
        pad_s_req=core_s-250; pad_e_req=core_e+250
        pad_complete=(pad_s_req>=0 and pad_e_req<=len(comb))
        pad_s=max(0,pad_s_req); pad_e=min(len(comb),pad_e_req)
        pad_seq=comb[pad_s:pad_e]
        core_hash=hashlib.sha256(core_seq.encode()).hexdigest(); pad_hash=hashlib.sha256(pad_seq.encode()).hexdigest()
        pseq=partner_fa.get(a['nearest_downstream_orf_id'],'')
        if not pseq: raise SystemExit(f'FAIL partner peptide missing {a["nearest_downstream_orf_id"]}')
        rt_hash=hashlib.sha256(pep.encode()).hexdigest(); partner_hash=hashlib.sha256(pseq.encode()).hexdigest()
        row={
          'family_rep_id':fam,'candidate_id':cid,'accession':l['accession'],'parent_contig_key':l['parent_contig_key'],
          'rt_strand':strand,'rt_aa_length':len(pep),'rt_peptide_sha256':rt_hash,
          'upstream_available_bp':l['upstream_available_bp'],'downstream_available_bp':l['downstream_available_bp'],
          'strict_array_n':len(arr_by[cid]),'strict_leftmost_repeat_start_rel0':left_rel,
          'dominant_recurrent_partner_cluster_id':dominant_cluster[fam],
          'partner_orf_id':a['nearest_downstream_orf_id'],'partner_aa_length':len(pseq),'partner_peptide_sha256':partner_hash,
          'partner_orientation':a['nearest_orf_orientation_relative_to_rt'],'partner_gap_bp':a['nearest_orf_gap_after_rt_stop_bp'],
          'architecture_core_length_bp':len(core_seq),'architecture_core_sha256':core_hash,
          'architecture_padded_250bp_complete':pad_complete,'architecture_padded_length_bp':len(pad_seq),'architecture_padded_sha256':pad_hash,
        }
        context_rows.append(row)
        core_fa.append(f'>{cid} family={fam} accession={l["accession"]} core=strict_repeat_to_partner_end\n{wrap(core_seq)}')
        padded_fa.append(f'>{cid} family={fam} accession={l["accession"]} padded250_complete={str(pad_complete).lower()}\n{wrap(pad_seq)}')
        rt_core_fa.append(f'>{cid} family={fam} accession={l["accession"]}\n{wrap(pep)}')
        partner_core_fa.append(f'>{a["nearest_downstream_orf_id"]} family={fam} candidate={cid} accession={l["accession"]}\n{wrap(pseq)}')

    # architecture-linked representative: max min flank, total flank, then candidate id lexical.
    reps=[]
    for fam in [r['family_rep_id'] for r in targets]:
        xs=[r for r in context_rows if r['family_rep_id']==fam]
        best=sorted(xs,key=lambda r:(-min(int(r['upstream_available_bp']),int(r['downstream_available_bp'])),-(int(r['upstream_available_bp'])+int(r['downstream_available_bp'])),r['candidate_id']))[0]
        reps.append(best)

    # Prior structure family coverage from ER1.
    er1=list(csv.DictReader(open(er1_sum,encoding='utf-8'),delimiter='\t'))
    existing={r['family_rep_id']:r for r in er1 if r['family_rep_id'] in famset}
    if len(existing)!=1:
        raise SystemExit(f'FAIL expected exactly one pre-existing structure-tested target family, got {len(existing)} {sorted(existing)}')
    existing_fam=next(iter(existing))
    newreps=[r for r in reps if r['family_rep_id'] not in existing]
    if len(newreps)!=6: raise SystemExit(f'FAIL new structure target n={len(newreps)}')

    # Write manifests / sequences.
    tf=out/'F3B4D_TARGET_FAMILIES_7.tsv'
    tcols=['priority_rank','family_rep_id','formal_locus_n','accession_n','strict_array_locus_n','strict_array_dominant_partner_cluster_id','strict_array_dominant_partner_cluster_instance_n','strict_array_dominant_partner_cluster_accession_n','strict_array_dominant_partner_cluster_median_gap_bp','preexisting_structure_tested','preexisting_structure_candidate_id','preexisting_structure_tier']
    with tf.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=tcols,delimiter='\t',lineterminator='\n'); w.writeheader()
        for i,t in enumerate(targets,1):
            ex=existing.get(t['family_rep_id'],{})
            w.writerow({
              'priority_rank':i,'family_rep_id':t['family_rep_id'],'formal_locus_n':t['formal_locus_n'],'accession_n':t['accession_n'],
              'strict_array_locus_n':t['strict_array_locus_n'],'strict_array_dominant_partner_cluster_id':t['strict_array_dominant_partner_cluster_id'],
              'strict_array_dominant_partner_cluster_instance_n':t['strict_array_dominant_partner_cluster_instance_n'],
              'strict_array_dominant_partner_cluster_accession_n':t['strict_array_dominant_partner_cluster_accession_n'],
              'strict_array_dominant_partner_cluster_median_gap_bp':t['strict_array_dominant_partner_cluster_median_gap_bp'],
              'preexisting_structure_tested':bool(ex),'preexisting_structure_candidate_id':ex.get('candidate_id',''),
              'preexisting_structure_tier':ex.get('b3er1_final_best_tier','') or ex.get('f3b3e_final_best_tier','') or ex.get('final_best_tier','')
            })
    cr=out/'F3B4D_CORE_RECURRENT_STRICT_LOCI_16.tsv'; ccols=list(context_rows[0])
    with cr.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=ccols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(sorted(context_rows,key=lambda r:(r['family_rep_id'],r['accession'],r['candidate_id'])))
    (inp/'F3B4D_CORE_ARCHITECTURE_CONTEXT_16.fna').write_text(''.join(core_fa),encoding='ascii')
    (inp/'F3B4D_PADDED250_ARCHITECTURE_CONTEXT_16.fna').write_text(''.join(padded_fa),encoding='ascii')
    (inp/'F3B4D_CORE_RT_16.faa').write_text(''.join(rt_core_fa),encoding='ascii')
    (inp/'F3B4D_CORE_PARTNER_16.faa').write_text(''.join(partner_core_fa),encoding='ascii')

    # 7 architecture representatives and 7 partner-cluster representatives.
    rt7=[]; partner7=[]; rep_rows=[]
    for r in reps:
        cid=r['candidate_id']; fam=r['family_rep_id']; p=im[cid]; pep=peps[p['orfdb_key']]
        pseq=partner_fa[r['partner_orf_id']]
        rt7.append(f'>{cid} family={fam} role=architecture_linked_representative\n{wrap(pep)}')
        partner7.append(f'>{r["partner_orf_id"]} family={fam} role=recurrent_strict_partner_cluster_member\n{wrap(pseq)}')
        z=dict(r); z['architecture_representative_rule']='max min(upstream,downstream), then total available context, then candidate_id lexical'; rep_rows.append(z)
    (inp/'F3B4D_RT_ARCHITECTURE_REPRESENTATIVES_7.faa').write_text(''.join(rt7),encoding='ascii')
    (inp/'F3B4D_PARTNER_REPRESENTATIVES_7.faa').write_text(''.join(partner7),encoding='ascii')
    rp=out/'F3B4D_ARCHITECTURE_REPRESENTATIVE_PER_FAMILY_7.tsv'; rcols=list(rep_rows[0])
    with rp.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=rcols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rep_rows)

    # New 6 structure targets + Boltz YAML.
    dst=out/'F3B4D_D_STRUCTURE_EXTENSION_TARGETS_6.tsv'
    dcols=['extension_rank','family_rep_id','candidate_id','accession','aa_length','yaml_path','sequence_sha256','selection_rule']
    dr=[]; dfa=[]
    for i,r in enumerate(newreps,1):
        cid=r['candidate_id']; p=im[cid]; pep=peps[p['orfdb_key']]
        y=ydir/f'{cid}.yaml'; y.write_text('version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: '+pep+'\n')
        dr.append({'extension_rank':i,'family_rep_id':r['family_rep_id'],'candidate_id':cid,'accession':r['accession'],'aa_length':len(pep),'yaml_path':str(y),'sequence_sha256':hashlib.sha256(pep.encode()).hexdigest(),'selection_rule':'strict+recurrent-partner locus; max min(upstream,downstream), then total context, then candidate_id lexical'})
        dfa.append(f'>{cid} family={r["family_rep_id"]} accession={r["accession"]}\n{wrap(pep)}')
    with dst.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=dcols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(dr)
    (inp/'F3B4D_D_STRUCTURE_EXTENSION_6.faa').write_text(''.join(dfa),encoding='ascii')

    # Freeze lock.
    lock={
      'stage':'Stage05F3B4D-implementation-freeze','gate':'PASS_F3B4D_SEVEN_FAMILY_A_TO_E_RULES_FROZEN_PRE_CLOSEOUT_RESULTS','gate_pass':True,
      'target_family_rule':'all and only F3B4C2 families with strict_array_plus_recurrent_partner_cluster=True','target_family_n':7,'target_family_ids':[r['family_rep_id'] for r in targets],
      'core_architecture_locus_rule':'within each target family, strict ART-geometry locus whose nearest downstream ORF belongs to that family strict-array dominant recurrent partner cluster','core_architecture_locus_n':16,
      'A_lineage_rule':'exact architecture-core SHA256 plus near-context graph using 31-mer Jaccard >=0.95 and length ratio >=0.95 within family; this is near-context dereplication, not evolutionary lineage proof','A_architecture_core':'leftmost strict repeat start through nearest partner ORF stop; 250-bp padding audited separately',
      'B_repeat_rule':'all strict arrays from the seven frozen families; motif-core neighborhoods use Levenshtein distance <=1; ART detector thresholds are not changed','B_repeat_detector_changed':False,
      'C_partner_rule':'nearest downstream ORF cluster already frozen by F3B4C; conservation evaluated only within the recurrent strict-array dominant cluster; nearest ORF remains a partner candidate, not confirmed dedicated partner','C_partner_catalog_role':'descriptive Pfam/PHROG/VOG annotation if local resources are present',
      'D_structure_rule':'all and only target families not previously structure-tested in F3B3E-R1; one architecture-linked representative per family chosen by technical context completeness only','D_new_structure_family_n':6,'D_existing_structure_family_id':existing_fam,
      'D_boltz_model':'boltz2','D_msa_mode':'server-assisted','D_seed':20261007,'D_recycling_steps':3,'D_sampling_steps':200,'D_diffusion_samples':1,'D_no_kernels':True,'D_serial_exclusive_idle_gpu':True,
      'D_r5c_strong_rule':'alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50','D_r5c_suggestive_rule':'alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30',
      'E_known_system_rule':'same local frozen known-system reference-panel construction and MMseqs settings as F3B3B; optional user-supplied global FASTA and full ART catalogue are handled separately','E_global_novelty_decision_authorized_without_global_reference':False,'E_anthropic_unreported_claim_authorized_without_full_ART95_catalogue':False,
      'final_priority_rule':'computational evidence board only: R5C STRONG>SUGGESTIVE>WEAK/NO_HIT>technical unresolved; then recurrent strict-partner accession count; then near-context cluster count; no novelty claim unless E explicitly closes it',
      'b4c2_gate_sha256':sha256_file(gate_p),'b4c_family_summary_sha256':sha256_file(famsum_p),'b4c_locus_architecture_sha256':sha256_file(arch_p),'b4c_arrays_sha256':sha256_file(arr_p),
      'r4_instance_map_sha256':sha256_file(imap_p),'er1_summary_sha256':sha256_file(er1_sum),'vendor_context_helper_sha256':sha256_file(vendor),
      'target_family_manifest_sha256':sha256_file(tf),'core_locus_manifest_sha256':sha256_file(cr),'structure_target_manifest_sha256':sha256_file(dst),
      'claim_boundary':'Seven-family post-discovery computational closeout. A/B/C/D/E strengthen architecture, structure and novelty auditing but do not establish biochemical RT activity, dedicated partner function, or Anthropic-unreported novelty by themselves.'
    }
    Path(lock_json).write_text(json.dumps(lock,indent=2)+'\n')
    print(json.dumps(lock,indent=2))
    print('\n=== TARGET FAMILIES 7 ===')
    with open(tf) as f: print(f.read().rstrip())
    print('\n=== STRUCTURE EXTENSION TARGETS 6 ===')
    with open(dst) as f: print(f.read().rstrip())
    print('PASS_F3B4D_PREPARE_AND_FREEZE_COMPLETE')

if __name__=='__main__': main()
