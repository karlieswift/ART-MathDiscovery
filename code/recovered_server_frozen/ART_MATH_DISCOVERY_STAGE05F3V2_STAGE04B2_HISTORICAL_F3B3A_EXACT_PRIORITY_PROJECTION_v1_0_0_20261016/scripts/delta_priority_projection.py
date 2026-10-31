#!/usr/bin/env python3
import argparse,csv,gzip,hashlib,json,math,re,sys
from collections import defaultdict,Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from rule_core import repkey,select_families,qualifying_profile_pass

MOTIF_RE=re.compile(r'(?=([YF][A-Z]DD))')
EXPECTED={
 'stage04a_gate_sha256':'5edb1445ed08c1120e855ec81dad10837d5e5b0a80b01d5def727716c4a8c88a',
 'stage04b1_gate_sha256':'a71a020736c9a7ac22db0102e3469228dda7a74c9555e6da4c6d05de8b74836e',
 'selected_manifest_sha256':'8faa369c8c974e01a3caccc6ed666210735472fcbd3b144714a662576d939eb5',
 'formal_instances_sha256':'c73d1852f7c7d3d82c6ad62d10872a3d1fd05fd67b477bf93706994faef03e97',
 'combined_family_membership_sha256':'fb84dca6ff346b2438b736b95a7dd606a7163e476e8d0e162ff70a72a6b3a0c0',
}

def shafile(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def shaseq(s): return hashlib.sha256(s.encode('ascii')).hexdigest()

def fasta_iter(p):
    h=None; seq=[]
    with open(p,encoding='ascii') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: yield h,''.join(seq).upper()
                h=line[1:].split()[0]; seq=[]
            else: seq.append(line)
        if h is not None: yield h,''.join(seq).upper()

def read_tsv(p):
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows,fields):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage04a-root',required=True)
    ap.add_argument('--stage04b1-root',required=True)
    ap.add_argument('--delta-root',required=True)
    ap.add_argument('--semantics-audit',required=True)
    ap.add_argument('--historical-replay-gate',required=True)
    ap.add_argument('--out',required=True)
    a=ap.parse_args()
    A=Path(a.stage04a_root).resolve(); B=Path(a.stage04b1_root).resolve(); D=Path(a.delta_root).resolve(); O=Path(a.out).resolve(); O.mkdir(parents=True,exist_ok=True)

    satp=A/'outputs/STAGE04A_60G_SATURATION_CENSUS.json'
    famp=A/'outputs/COMBINED_FAMILY_SUMMARY.tsv'
    memp=A/'outputs/COMBINED_FAMILY_MEMBERSHIP.tsv'
    mapglob=A/'outputs/GLOBAL_EXACT_AA_MAP.tsv'
    b1p=B/'outputs/STAGE04B1_DELTA_FAMILY_REPLICATION_CENSUS.json'
    b1fam=B/'outputs/DELTA_REPRESENTED_FAMILY_REPLICATION_CENSUS.tsv'
    manp=D/'frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'
    finst=D/'outputs/F3V2_FORMAL_INSTANCES.tsv'
    dg=D/'outputs/STAGE05F3V2_GATE.json'
    sem=Path(a.semantics_audit); hrep=Path(a.historical_replay_gate)
    for p in [satp,famp,memp,mapglob,b1p,b1fam,manp,finst,dg,sem,hrep]:
        if not p.is_file(): raise RuntimeError('Missing required input: '+str(p))

    if shafile(satp)!=EXPECTED['stage04a_gate_sha256']: raise RuntimeError('Stage04A gate SHA changed')
    if shafile(b1p)!=EXPECTED['stage04b1_gate_sha256']: raise RuntimeError('Stage04B1 gate SHA changed')
    if shafile(manp)!=EXPECTED['selected_manifest_sha256']: raise RuntimeError('Delta selected manifest SHA changed')
    if shafile(finst)!=EXPECTED['formal_instances_sha256']: raise RuntimeError('Delta formal instances SHA changed')
    if shafile(memp)!=EXPECTED['combined_family_membership_sha256']: raise RuntimeError('Combined family membership SHA changed')

    sat=json.load(open(satp)); b1=json.load(open(b1p)); delta_gate=json.load(open(dg)); sa=json.load(open(sem)); hr=json.load(open(hrep))
    if sat.get('gate_pass') is not True or b1.get('gate_pass') is not True or delta_gate.get('gate_pass') is not True: raise RuntimeError('Upstream gate not PASS')
    if sa.get('gate_pass') is not True: raise RuntimeError('Historical model-semantics audit not PASS')
    if hr.get('gate_pass') is not True: raise RuntimeError('Historical F3B3A replay not PASS')
    qdef=sa.get('qualifying_definition_used_for_delta')
    if not qdef: raise RuntimeError('No validated qualifying definition')

    famsummary=read_tsv(famp)
    famclass={r['family_rep']:r['family_class'] for r in famsummary}
    delta_only_fams={fr for fr,c in famclass.items() if c=='DELTA_ONLY'}
    if len(delta_only_fams)!=24810: raise RuntimeError(f'DELTA_ONLY family n changed: {len(delta_only_fams)}')

    sha_to_family={}; sha_to_member={}
    with open(memp,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            s=r['sequence_sha256']; fr=r['family_rep']
            if s in sha_to_family and sha_to_family[s]!=fr: raise RuntimeError('Exact SHA maps to multiple families: '+s)
            sha_to_family[s]=fr; sha_to_member[s]=r['member']

    b1rows=read_tsv(b1fam)
    b1delta={r['family_rep']:r for r in b1rows if r['family_class']=='DELTA_ONLY'}
    if set(b1delta)!=delta_only_fams: raise RuntimeError('Stage04B1 DELTA_ONLY family set mismatch')

    # Frozen formal instance rows, grouped by accession.
    byacc=defaultdict(dict); formal_total=0
    with open(finst,newline='') as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            if row.get('formal')!='True': continue
            acc=row['source_accession']; q=row['query_id']
            if q in byacc[acc]: raise RuntimeError('Duplicate formal query within accession: '+q)
            byacc[acc][q]=row; formal_total+=1
    if formal_total!=88487: raise RuntimeError(f'Formal instance n changed: {formal_total}')

    mans=read_tsv(manp); accs=[r['accession'] for r in mans]
    if len(accs)!=254 or len(set(accs))!=254: raise RuntimeError('Delta accession universe changed')
    if set(byacc)!=set(accs): raise RuntimeError(f'Formal accession set mismatch formal={len(byacc)} manifest={len(accs)}')

    cand=[]; byfam=defaultdict(list); observed_formal=0; excluded_mixed=0
    for ai,acc in enumerate(accs,1):
        formal_rows=byacc[acc]
        faa=D/'outputs/formal_faa'/f'{acc}.formal.faa'
        dom=D/'work/per_accession'/acc/'RVT_ALL.domtblout'
        if not faa.is_file() or not dom.is_file(): raise RuntimeError('Missing per-accession frozen output for '+acc)
        seq={q:s for q,s in fasta_iter(faa)}
        if set(seq)!=set(formal_rows):
            raise RuntimeError(f'Formal FASTA/query mismatch for {acc}: fasta={len(seq)} table={len(formal_rows)}')

        formal_prof=defaultdict(set); qual_prof=defaultdict(set); bestdom={}
        with open(dom,errors='ignore') as f:
            for line in f:
                if not line.strip() or line.startswith('#'): continue
                x=line.split()
                if len(x)<22: continue
                profile=x[0]; q=x[3]
                if q not in formal_rows: continue
                full=float(x[6]); ie=float(x[12]); score=float(x[13])
                ali_from=int(x[17]); ali_to=int(x[18]); env_from=int(x[19]); env_to=int(x[20])
                if full<=1e-4 and ie<=1e-4:
                    formal_prof[q].add(profile)
                    rec={'profile':profile,'full_evalue':full,'domain_i_evalue':ie,'score':score,'ali_from':ali_from,'ali_to':ali_to,'env_from':env_from,'env_to':env_to}
                    key=(ie,full,profile,ali_from,ali_to)
                    old=bestdom.get(q)
                    if old is None or key < (old['domain_i_evalue'],old['full_evalue'],old['profile'],old['ali_from'],old['ali_to']): bestdom[q]=rec
                if qualifying_profile_pass(qdef,full,ie): qual_prof[q].add(profile)

        for q,row in formal_rows.items():
            observed_formal+=1
            s=seq[q]; ss=shaseq(s); fr=sha_to_family.get(ss)
            if fr is None: raise RuntimeError('Formal AA SHA absent from combined family map: '+ss)
            if fr not in delta_only_fams:
                excluded_mixed+=1; continue
            d=bestdom.get(q)
            if d is None: raise RuntimeError('Formal query lacks formal domtblout row: '+q)
            # Frozen scan best i-E must agree with direct formal-domtblout reconstruction.
            if float(row['best_domain_i_evalue']) != float(d['domain_i_evalue']):
                raise RuntimeError(f'Best domain i-E mismatch {q}: scan={row["best_domain_i_evalue"]} dom={d["domain_i_evalue"]}')
            lo=max(1,d['env_from']); hi=min(len(s),d['env_to']); sub=s[lo-1:hi]
            hits=[(lo+m.start(),m.group(1)) for m in MOTIF_RE.finditer(sub)]
            ghits=[(1+m.start(),m.group(1)) for m in MOTIF_RE.finditer(s)]
            c={
              'candidate_id':q,'original_query_name':q,'family_rep_id':fr,'exact_rep_id':sha_to_member[ss],
              'accession':acc,'aa_length':len(s),'formal_model_n':len(formal_prof[q]),'qualifying_model_n':len(qual_prof[q]),
              'best_full_evalue':float(row['full_evalue']),'best_domain_i_evalue':float(row['best_domain_i_evalue']),
              'best_formal_profile':d['profile'],'best_domain_env_from':lo,'best_domain_env_to':hi,
              'motif_c_like_in_best_domain':bool(hits),'motif_c_like_domain_n':len(hits),
              'motif_c_like_domain_first_pos':hits[0][0] if hits else '',
              'motif_c_like_domain_first_seq':hits[0][1] if hits else '',
              'motif_c_like_anywhere':bool(ghits),'motif_c_like_anywhere_n':len(ghits),'sequence_sha256':ss,
            }
            cand.append(c); byfam[fr].append(c)
        if ai%25==0 or ai==len(accs): print(f'[DELTA PRIORITY] accession {ai}/{len(accs)} candidates_kept={len(cand)}',flush=True)

    if observed_formal!=88487: raise RuntimeError('Observed formal count mismatch')
    if set(byfam)!=delta_only_fams: raise RuntimeError(f'DELTA_ONLY family coverage mismatch missing={len(delta_only_fams-set(byfam))} extra={len(set(byfam)-delta_only_fams)}')

    # Family-level integrity against Stage04B1 descriptors.
    famrows=[]; family_integrity_bad=[]
    for fr in sorted(byfam):
        members=byfam[fr]; rep=sorted(members,key=repkey)[0]
        accn=len({x['accession'] for x in members}); exactn=len({x['sequence_sha256'] for x in members})
        b=b1delta[fr]
        observed=(len(members),accn,exactn)
        expected=(int(b['n_formal_instances']),int(b['n_accessions']),int(b['n_delta_exact_sequences']))
        if observed!=expected:
            if len(family_integrity_bad)<20: family_integrity_bad.append((fr,expected,observed))
            continue
        motif_n=sum(bool(x['motif_c_like_in_best_domain']) for x in members)
        famrows.append({
          'family_rep_id':fr,'candidate_instance_n':len(members),'accession_n':accn,'multirun_replication':accn>=2,
          'max_aa_length':max(x['aa_length'] for x in members),'max_formal_model_n':max(x['formal_model_n'] for x in members),
          'max_qualifying_model_n':max(x['qualifying_model_n'] for x in members),'motif_c_like_candidate_n':motif_n,
          'motif_c_like_family':motif_n>0,'priority_candidate_id':rep['candidate_id'],'priority_candidate_accession':rep['accession'],
          'priority_candidate_aa_length':rep['aa_length'],'priority_candidate_formal_model_n':rep['formal_model_n'],
          'priority_candidate_best_domain_i_evalue':rep['best_domain_i_evalue'],'priority_candidate_motif_c_like':rep['motif_c_like_in_best_domain'],
          'priority_candidate_best_formal_profile':rep['best_formal_profile'],
        })
    if family_integrity_bad or len(famrows)!=24810:
        raise RuntimeError(f'Family metric integrity failure bad_n={len(family_integrity_bad)} examples={family_integrity_bad[:5]} famrows={len(famrows)}')

    selected,reason,rank,fam_by_id=select_families(famrows)
    for x in famrows:
        fr=x['family_rep_id']; x['shortlisted']=fr in rank; x['shortlist_rank']=rank.get(fr,''); x['shortlist_reason']=';'.join(reason.get(fr,[]))
    famrows.sort(key=lambda x:(0 if x['shortlisted'] else 1, x['shortlist_rank'] if x['shortlisted'] else 10**9, -x['accession_n'],-x['max_aa_length'],x['family_rep_id']))

    cand_fields=list(cand[0]); fam_fields=list(famrows[0])
    write_tsv(O/'F3V2_STAGE04B2_CANDIDATE_CHARACTERIZATION_DELTA_ONLY.tsv',cand,cand_fields)
    write_tsv(O/'F3V2_STAGE04B2_FAMILY_PRIORITY_DELTA_ONLY.tsv',famrows,fam_fields)
    shortlist=[]
    for fr in selected:
        x=dict(fam_by_id[fr]); x.update(shortlisted=True,shortlist_rank=rank[fr],shortlist_reason=';'.join(reason[fr])); shortlist.append(x)
    write_tsv(O/'F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.tsv',shortlist,fam_fields)

    seq_by_q={x['candidate_id']:None for x in cand if x['candidate_id'] in {fam_by_id[fr]['priority_candidate_id'] for fr in selected}}
    # Avoid retaining all sequences: re-read only selected accessions/queries.
    need=defaultdict(set)
    for fr in selected:
        x=fam_by_id[fr]; need[x['priority_candidate_accession']].add(x['priority_candidate_id'])
    seqsel={}
    for acc,qs in need.items():
        for q,s in fasta_iter(D/'outputs/formal_faa'/f'{acc}.formal.faa'):
            if q in qs: seqsel[q]=s
    wantq={fam_by_id[fr]['priority_candidate_id'] for fr in selected}
    if set(seqsel)!=wantq: raise RuntimeError(f'Shortlist FASTA recovery mismatch missing={wantq-set(seqsel)}')
    with open(O/'F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.faa','w') as f:
        for fr in selected:
            x=fam_by_id[fr]; q=x['priority_candidate_id']; s=seqsel[q]
            f.write(f'>{q} family={fr} accession={x["priority_candidate_accession"]} aa_length={len(s)}\n')
            for i in range(0,len(s),80): f.write(s[i:i+80]+'\n')

    axis_counts=Counter()
    for fr in selected:
        for rr in reason[fr]: axis_counts[rr]+=1
    gate={
      'stage':'Stage05F3V2-Stage04B2-historical-F3B3A-exact-priority-projection',
      'gate':'PASS_F3V2_STAGE04B2_HISTORICAL_F3B3A_EXACT_RULE_PROJECTED_TO_DELTA_ONLY_FAMILIES',
      'gate_pass':True,
      'delta_formal_instance_n':88487,'delta_only_formal_instance_n':len(cand),'mixed_family_formal_instance_excluded_n':excluded_mixed,
      'delta_only_combined_family_n':len(famrows),'motif_c_like_candidate_n':sum(bool(x['motif_c_like_in_best_domain']) for x in cand),
      'motif_c_like_family_n':sum(bool(x['motif_c_like_family']) for x in famrows),'structure_priority_shortlist_family_n':len(selected),
      'historical_rule_source_stage':'Stage05F3B3A (2026-10-11)','historical_rule_replay_pass':True,
      'historical_model_semantics_reconstruction_pass':True,'qualifying_model_definition':qdef,
      'shortlist_rule':'exact historical four-axis fixed quota (8 replication, 8 max length, 8 model coverage, 8 motif topology) with deterministic backfill to 32 unique families',
      'shortlist_axis_reason_counts':dict(axis_counts),'family_universe':'Stage04A analysis-only DELTA_ONLY combined families',
      'known_rt_annotation_used_for_ranking':False,'novelty_used_for_ranking':False,'structure_used_for_ranking':False,'truth_used_for_ranking':False,
      'retuning_authorized':False,'scientific_threshold_changed':False,'science_change':False,
      'candidate_characterization_sha256':shafile(O/'F3V2_STAGE04B2_CANDIDATE_CHARACTERIZATION_DELTA_ONLY.tsv'),
      'family_priority_sha256':shafile(O/'F3V2_STAGE04B2_FAMILY_PRIORITY_DELTA_ONLY.tsv'),
      'shortlist_tsv_sha256':shafile(O/'F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.tsv'),
      'shortlist_fasta_sha256':shafile(O/'F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.faa'),
      'claim_boundary':'The ranking logic is the exact historical F3B3A compute-priority rule frozen before the 54-Gb delta existed. It is projected onto the Stage04A analysis-only DELTA_ONLY combined-family universe. This is compute triage, not evidence of RT biochemical activity, ART identity, biological independence, or novelty.'
    }
    with open(O/'STAGE04B2_GATE.json','w') as f: json.dump(gate,f,indent=2); f.write('\n')
    outs=[O/'F3V2_STAGE04B2_CANDIDATE_CHARACTERIZATION_DELTA_ONLY.tsv',O/'F3V2_STAGE04B2_FAMILY_PRIORITY_DELTA_ONLY.tsv',O/'F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.tsv',O/'F3V2_STAGE04B2_STRUCTURE_PRIORITY_SHORTLIST_32.faa',O/'STAGE04B2_GATE.json']
    with open(O/'SHA256SUMS.txt','w') as f:
        for p in outs: f.write(f'{shafile(p)}  {p.name}\n')
    print(json.dumps(gate,indent=2))
    print('\n=== DELTA STRUCTURE PRIORITY SHORTLIST 32 ===')
    print('\t'.join(['rank','family_rep_id','reason','accession_n','candidate_instance_n','max_aa_length','max_formal_model_n','motif_n','priority_candidate_id','accession']))
    for fr in selected:
        x=fam_by_id[fr]
        print('\t'.join(map(str,[rank[fr],fr,';'.join(reason[fr]),x['accession_n'],x['candidate_instance_n'],x['max_aa_length'],x['max_formal_model_n'],x['motif_c_like_candidate_n'],x['priority_candidate_id'],x['priority_candidate_accession']])))
    print('PASS_STAGE05F3V2_STAGE04B2_HISTORICAL_F3B3A_EXACT_PRIORITY_PROJECTION_COMPLETE')

if __name__=='__main__':
    try: main()
    except Exception as e:
        print('[HOLD]',e,file=sys.stderr)
        raise
