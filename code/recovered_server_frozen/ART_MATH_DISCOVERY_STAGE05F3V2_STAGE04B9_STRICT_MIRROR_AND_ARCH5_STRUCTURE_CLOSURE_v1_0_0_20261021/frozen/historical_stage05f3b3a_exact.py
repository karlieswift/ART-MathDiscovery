#!/usr/bin/env python3
import csv, gzip, hashlib, json, re, sys
from collections import defaultdict
from pathlib import Path


def sha256_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''):
            h.update(b)
    return h.hexdigest()


def read_fasta(path):
    h=None; seq=[]
    with open(path,encoding='ascii') as f:
        for line in f:
            line=line.rstrip('\r\n')
            if line.startswith('>'):
                if h is not None:
                    yield h,''.join(seq).upper()
                h=line[1:]; seq=[]
            elif line:
                seq.append(line.strip())
        if h is not None:
            yield h,''.join(seq).upper()


def ffloat(x):
    try: return float(x)
    except: return float('inf')


def main():
    if len(sys.argv)!=11:
        raise SystemExit('usage: stage05f3b3a.py R4_GATE R4_CANDIDATE_FAMILY_MAP R4_FAMILY_SUMMARY R3_FORMAL_IDS R3_NAMESPACE_FASTA R3_PAIR_HITS_GZ R3_DOMTBLOUT LOCK OUTDIR SHORTLIST_N')
    (r4_gate_file,fmap_file,famsum_file,formal_ids_file,namespace_fa,pair_gz,dom_file,lock_file,outdir,shortlist_n_s)=sys.argv[1:11]
    shortlist_n=int(shortlist_n_s)
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    lock=json.load(open(lock_file)); g=json.load(open(r4_gate_file))
    assert g['gate_pass'] is True
    assert g['formal_candidate_instance_n']==10598
    assert g['exact_aa_unique_n']==10329
    assert g['sequence_family_n']==4912
    assert g['selected_namespace_unknown_accession_n']==0
    assert g['formal_unknown_accession_n']==0
    assert g['family_membership_changed'] is False
    assert g['hmm_results_changed'] is False
    assert shortlist_n==lock['structure_priority_shortlist_family_n']==32

    formal=[x.strip() for x in open(formal_ids_file) if x.strip()]
    fset=set(formal)
    if len(formal)!=10598 or len(fset)!=10598:
        raise SystemExit(f'FAIL formal IDs n={len(formal)} unique={len(fset)}')

    # Frozen R4 family map with canonical accessions.
    fmap={}
    with open(fmap_file,encoding='utf-8') as f:
        r=csv.DictReader(f,delimiter='\t')
        req={'candidate_id','original_query_name','family_rep_id','exact_rep_id','aa_length','accession','formal_model_n','qualifying_model_n','best_full_evalue','best_domain_i_evalue'}
        if not req.issubset(set(r.fieldnames or [])):
            raise SystemExit(f'FAIL family map missing {sorted(req-set(r.fieldnames or []))}')
        for row in r:
            fmap[row['candidate_id']]=row
    if set(fmap)!=fset:
        raise SystemExit(f'FAIL R4 family map candidate set mismatch map={len(fmap)} formal={len(fset)}')

    famsum={}
    with open(famsum_file,encoding='utf-8') as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            famsum[row['family_rep_id']]=row
    if len(famsum)!=4912:
        raise SystemExit(f'FAIL family summary n={len(famsum)}')

    # Recover exact formal sequences from the namespace FASTA.
    seq={}
    for h,s in read_fasta(namespace_fa):
        q=h.split()[0]
        if q in fset:
            if q in seq: raise SystemExit(f'FAIL duplicate formal instance in namespace FASTA: {q}')
            seq[q]=s
    if set(seq)!=fset:
        raise SystemExit(f'FAIL formal sequence recovery missing={len(fset-set(seq))}')

    # Profile support summary; descriptive only.
    prof=defaultdict(list)
    with gzip.open(pair_gz,'rt') as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            q=row['query_id']
            if q in fset:
                prof[q].append(row)
    for q in formal:
        if not prof[q] or not any(x['class']=='FORMAL' for x in prof[q]):
            raise SystemExit(f'FAIL formal candidate lacks formal profile evidence: {q}')

    # Best formal domain from the unchanged R3 hmmscan. HMMER domtblout columns:
    # target profile ... query ... full-E ... i-E ... ali-from ali-to env-from env-to ...
    bestdom={}
    with open(dom_file,encoding='utf-8') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            x=line.split()
            if len(x)<22: continue
            profile=x[0]; q=x[3]
            if q not in fset: continue
            full_e=float(x[6]); ie=float(x[12])
            if full_e>1e-4 or ie>1e-4: continue
            rec={
                'profile':profile,'full_evalue':full_e,'domain_i_evalue':ie,
                'ali_from':int(x[17]),'ali_to':int(x[18]),
                'env_from':int(x[19]),'env_to':int(x[20])
            }
            old=bestdom.get(q)
            key=(ie,full_e,profile,rec['ali_from'],rec['ali_to'])
            if old is None or key < (old['domain_i_evalue'],old['full_evalue'],old['profile'],old['ali_from'],old['ali_to']):
                bestdom[q]=rec
    missing=fset-set(bestdom)
    if missing:
        raise SystemExit(f'FAIL formal candidates without formal domtblout row n={len(missing)} examples={sorted(missing)[:20]}')

    motif_re=re.compile(r'(?=([YF][A-Z]DD))')
    cand=[]
    byfam=defaultdict(list)
    for q in formal:
        r=fmap[q]; s=seq[q]; d=bestdom[q]
        # 1-based HMMER coordinates; inspect best-domain envelope only.
        lo=max(1,d['env_from']); hi=min(len(s),d['env_to'])
        sub=s[lo-1:hi]
        hits=[]
        for m in motif_re.finditer(sub):
            pos=lo+m.start()
            hits.append((pos,m.group(1)))
        global_hits=[]
        for m in motif_re.finditer(s):
            global_hits.append((1+m.start(),m.group(1)))
        c={
            'candidate_id':q,
            'original_query_name':r['original_query_name'],
            'family_rep_id':r['family_rep_id'],
            'exact_rep_id':r['exact_rep_id'],
            'accession':r['accession'],
            'aa_length':int(r['aa_length']),
            'formal_model_n':int(r['formal_model_n']),
            'qualifying_model_n':int(r['qualifying_model_n']),
            'best_full_evalue':float(r['best_full_evalue']),
            'best_domain_i_evalue':float(r['best_domain_i_evalue']),
            'best_formal_profile':d['profile'],
            'best_domain_env_from':lo,
            'best_domain_env_to':hi,
            'motif_c_like_in_best_domain':bool(hits),
            'motif_c_like_domain_n':len(hits),
            'motif_c_like_domain_first_pos':hits[0][0] if hits else '',
            'motif_c_like_domain_first_seq':hits[0][1] if hits else '',
            'motif_c_like_anywhere':bool(global_hits),
            'motif_c_like_anywhere_n':len(global_hits),
            'sequence_sha256':hashlib.sha256(s.encode('ascii')).hexdigest(),
        }
        cand.append(c); byfam[c['family_rep_id']].append(c)

    ccols=list(cand[0])
    with open(out/'F3B3A_CANDIDATE_CHARACTERIZATION.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=ccols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(cand)

    # Representative candidate per family for expensive orthogonal follow-up.
    def repkey(c):
        return (
            0 if c['motif_c_like_in_best_domain'] else 1,
            -c['formal_model_n'],
            -c['aa_length'],
            c['best_domain_i_evalue'],
            c['candidate_id']
        )
    famrows=[]
    for fr, members in byfam.items():
        rep=sorted(members,key=repkey)[0]
        fs=famsum[fr]
        motif_n=sum(bool(x['motif_c_like_in_best_domain']) for x in members)
        famrows.append({
            'family_rep_id':fr,
            'candidate_instance_n':int(fs['candidate_instance_n']),
            'accession_n':int(fs['accession_n']),
            'multirun_replication':str(fs['multirun_replication']).lower()=='true',
            'max_aa_length':int(float(fs['max_aa_length'])),
            'max_formal_model_n':int(fs['max_formal_model_n']),
            'max_qualifying_model_n':int(fs['max_qualifying_model_n']),
            'motif_c_like_candidate_n':motif_n,
            'motif_c_like_family':motif_n>0,
            'priority_candidate_id':rep['candidate_id'],
            'priority_candidate_accession':rep['accession'],
            'priority_candidate_aa_length':rep['aa_length'],
            'priority_candidate_formal_model_n':rep['formal_model_n'],
            'priority_candidate_best_domain_i_evalue':rep['best_domain_i_evalue'],
            'priority_candidate_motif_c_like':rep['motif_c_like_in_best_domain'],
            'priority_candidate_best_formal_profile':rep['best_formal_profile'],
        })

    # Frozen 4-axis quota: 8 families each for replication, length, model coverage, motif topology.
    fam_by_id={x['family_rep_id']:x for x in famrows}
    quota=lock['shortlist_quota_per_axis']
    assert quota==8
    axis_order=lock['shortlist_axis_order']
    selected=[]; reason=defaultdict(list)
    def add(items,label):
        for x in items[:quota]:
            fr=x['family_rep_id']; reason[fr].append(label)
            if fr not in selected: selected.append(fr)
    add(sorted(famrows,key=lambda x:(-x['accession_n'],-x['candidate_instance_n'],-x['max_aa_length'],x['family_rep_id'])),'REPLICATION_TOP8')
    add(sorted(famrows,key=lambda x:(-x['max_aa_length'],-x['accession_n'],-x['max_formal_model_n'],x['family_rep_id'])),'LENGTH_TOP8')
    add(sorted(famrows,key=lambda x:(-x['max_formal_model_n'],-x['max_qualifying_model_n'],-x['accession_n'],-x['max_aa_length'],x['family_rep_id'])),'MODEL_COVERAGE_TOP8')
    motif_sorted=[x for x in famrows if x['motif_c_like_family']]
    motif_sorted=sorted(motif_sorted,key=lambda x:(-x['motif_c_like_candidate_n'],x['priority_candidate_best_domain_i_evalue'],-x['priority_candidate_aa_length'],-x['accession_n'],x['family_rep_id']))
    add(motif_sorted,'MOTIF_TOPOLOGY_TOP8')

    # Deterministic backfill to exactly 32 unique families when quotas overlap or motif families <8.
    backfill=sorted(famrows,key=lambda x:(
        0 if x['motif_c_like_family'] else 1,
        -x['accession_n'],-x['max_formal_model_n'],-x['max_aa_length'],
        -x['candidate_instance_n'],x['priority_candidate_best_domain_i_evalue'],x['family_rep_id']))
    for x in backfill:
        if len(selected)>=shortlist_n: break
        fr=x['family_rep_id']
        if fr not in selected:
            selected.append(fr); reason[fr].append('DETERMINISTIC_BACKFILL')
    if len(selected)!=shortlist_n:
        raise SystemExit(f'FAIL shortlist family n={len(selected)} expected={shortlist_n}')

    selected_rank={fr:i+1 for i,fr in enumerate(selected)}
    for x in famrows:
        fr=x['family_rep_id']
        x['shortlisted']=fr in selected_rank
        x['shortlist_rank']=selected_rank.get(fr,'')
        x['shortlist_reason']=';'.join(reason.get(fr,[]))
    famrows.sort(key=lambda x:(0 if x['shortlisted'] else 1, x['shortlist_rank'] if x['shortlisted'] else 10**9, -x['accession_n'],-x['max_aa_length'],x['family_rep_id']))
    fcols=list(famrows[0])
    with open(out/'F3B3A_FAMILY_PRIORITY.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fcols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(famrows)
    with open(out/'F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fcols,delimiter='\t',lineterminator='\n');w.writeheader();
        for fr in selected:
            w.writerow(fam_by_id[fr] | {'shortlisted':True,'shortlist_rank':selected_rank[fr],'shortlist_reason':';'.join(reason[fr])})

    # FASTA of 32 priority representatives.
    with open(out/'F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.faa','w') as f:
        for fr in selected:
            c=fam_by_id[fr]; q=c['priority_candidate_id']; s=seq[q]
            f.write(f'>{q} family={fr} accession={c["priority_candidate_accession"]} aa_length={len(s)}\n')
            for i in range(0,len(s),80): f.write(s[i:i+80]+'\n')

    gate={
        'stage':'Stage05F3B3A',
        'gate':'PASS_F3B3A_POSTFREEZE_PRIORITY_AND_CATALYTIC_TOPOLOGY',
        'gate_pass':True,
        'formal_candidate_instance_n':len(cand),
        'sequence_family_n':len(famrows),
        'formal_accession_n':g['formal_accession_n'],
        'multirun_family_n':g['multirun_family_n'],
        'max_accession_replication_n':g['max_accession_replication_n'],
        'motif_definition':'descriptive canonical motif-C-like [YF]xDD; not a biochemical identity threshold',
        'motif_c_like_candidate_n':sum(bool(x['motif_c_like_in_best_domain']) for x in cand),
        'motif_c_like_family_n':sum(bool(x['motif_c_like_family']) for x in famrows),
        'structure_priority_shortlist_family_n':len(selected),
        'shortlist_rule':'four-axis fixed quota (8 replication, 8 max length, 8 model coverage, 8 motif topology) with deterministic backfill to 32 unique families',
        'shortlist_is_inferential_threshold':False,
        'known_rt_annotation_used_for_ranking':False,
        'novelty_used_for_ranking':False,
        'candidate_membership_changed':False,
        'family_membership_changed':False,
        'hmm_results_changed':False,
        'scientific_threshold_changed':False,
        'r4_gate_sha256':sha256_file(r4_gate_file),
        'candidate_characterization_sha256':sha256_file(out/'F3B3A_CANDIDATE_CHARACTERIZATION.tsv'),
        'family_priority_sha256':sha256_file(out/'F3B3A_FAMILY_PRIORITY.tsv'),
        'shortlist_tsv_sha256':sha256_file(out/'F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv'),
        'shortlist_fasta_sha256':sha256_file(out/'F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.faa'),
        'claim_boundary':'F3B3A is post-freeze descriptive characterization and compute-priority triage. Motif-C-like presence, replication, ORF length, and profile breadth are priority attributes only; they do not establish RT/ART identity or novelty.'
    }
    with open(out/'STAGE05F3B3A_GATE.json','w') as f:
        json.dump(gate,f,indent=2);f.write('\n')
    print(json.dumps(gate,indent=2))

if __name__=='__main__': main()
