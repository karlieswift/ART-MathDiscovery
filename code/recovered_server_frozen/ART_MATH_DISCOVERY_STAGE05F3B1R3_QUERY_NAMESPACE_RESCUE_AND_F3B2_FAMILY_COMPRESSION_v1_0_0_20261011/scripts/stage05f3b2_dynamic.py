#!/usr/bin/env python3
import csv,gzip,hashlib,json,statistics,sys
from collections import defaultdict
from pathlib import Path


def sha256_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()


def read_fasta(path):
    h=None;seq=[]
    with open(path,encoding='ascii') as f:
        for line in f:
            line=line.rstrip('\r\n')
            if line.startswith('>'):
                if h is not None: yield h,''.join(seq).upper()
                h=line[1:];seq=[]
            elif line: seq.append(line.strip())
        if h is not None: yield h,''.join(seq).upper()


def ordkey(q,meta): return int(meta[q]['serialized_record_ordinal'])


def prepare(args):
    formal_file,numfa,map_tsv,pairgz,outdir=args
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    formal=[x.strip() for x in open(formal_file) if x.strip()]
    fset=set(formal)
    if not formal or len(fset)!=len(formal): raise SystemExit('FAIL formal instance IDs empty/duplicate')

    meta={}
    with open(map_tsv) as f:
        for r in csv.DictReader(f,delimiter='\t'): meta[r['instance_id']]=r
    if not fset <= set(meta): raise SystemExit(f'FAIL formal IDs missing from namespace map n={len(fset-set(meta))}')

    records={}
    for h,s in read_fasta(numfa):
        q=h.split()[0]
        if q in records: raise SystemExit(f'FAIL duplicate FASTA instance id {q}')
        if q in fset: records[q]=s
    if set(records)!=fset: raise SystemExit(f'FAIL formal FASTA membership mismatch missing={len(fset-set(records))} extra={len(set(records)-fset)}')

    pairs=defaultdict(list)
    with gzip.open(pairgz,'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            q=r['query_id']
            if q not in fset: continue
            pairs[q].append({'profile':r['profile'],'full_evalue':float(r['full_evalue']),'domain_ie':float(r['best_domain_i_evalue']),'class':r['class']})
    for q in formal:
        if not pairs[q] or not any(x['class']=='FORMAL' for x in pairs[q]): raise SystemExit(f'FAIL formal instance lacks formal HMM evidence {q}')

    seq_groups=defaultdict(list)
    for q,s in records.items(): seq_groups[s].append(q)
    candidate_to_exact={};exact_groups={}
    for seq,members in seq_groups.items():
        members=sorted(members,key=lambda q:ordkey(q,meta));rep=members[0];exact_groups[rep]=members
        for q in members:candidate_to_exact[q]=rep
    exact_reps=sorted(exact_groups,key=lambda q:ordkey(q,meta))

    with open(out/'F3B2_EXACT_AA_MAP.tsv','w') as f:
        f.write('candidate_id\texact_rep_id\texact_group_size\tsequence_sha256\n')
        for rep in exact_reps:
            sh=hashlib.sha256(records[rep].encode()).hexdigest()
            for q in exact_groups[rep]: f.write(f'{q}\t{rep}\t{len(exact_groups[rep])}\t{sh}\n')
    with open(out/'F3B2_EXACT_AA_REPRESENTATIVES.faa','w') as f:
        for rep in exact_reps:
            s=records[rep];f.write(f'>{rep}\n')
            for i in range(0,len(s),80):f.write(s[i:i+80]+'\n')

    with open(out/'F3B2_CANDIDATE_BASE.tsv','w') as f:
        f.write('candidate_id\toriginal_query_name\tserialized_record_ordinal\torfdb_key\taa_length\taccession\texact_rep_id\texact_group_size\tformal_model_n\tqualifying_model_n\tbest_full_evalue\tbest_domain_i_evalue\tsource_header\n')
        for q in sorted(formal,key=lambda q:ordkey(q,meta)):
            pp=pairs[q];rep=candidate_to_exact[q];m=meta[q]
            fm={x['profile'] for x in pp if x['class']=='FORMAL'};qm={x['profile'] for x in pp}
            vals=[q,m['original_query_name'],m['serialized_record_ordinal'],m['orfdb_key'],len(records[q]),m['accession'],rep,len(exact_groups[rep]),len(fm),len(qm),min(x['full_evalue'] for x in pp),min(x['domain_ie'] for x in pp),m['source_header']]
            f.write('\t'.join(map(str,vals))+'\n')
    summary={'formal_candidate_instance_n':len(formal),'exact_aa_unique_n':len(exact_reps),'exact_duplicate_instance_n':len(formal)-len(exact_reps)}
    with open(out/'F3B2_PREPARE_SUMMARY.json','w') as f:json.dump(summary,f,indent=2);f.write('\n')
    print(json.dumps(summary,indent=2));print('PASS_F3B2_EXACT_AA_COMPRESSION')


def finalize(args):
    base_file,exact_file,family_file,lock_file,outdir=args
    out=Path(outdir);lock=json.load(open(lock_file))
    base={}
    with open(base_file) as f:
        for r in csv.DictReader(f,delimiter='\t'): base[r['candidate_id']]=r
    exact={}
    with open(exact_file) as f:
        for r in csv.DictReader(f,delimiter='\t'):exact[r['candidate_id']]=r['exact_rep_id']
    exact_reps=set(exact.values());exact_to_family={}
    with open(family_file) as f:
        for line in f:
            if not line.strip():continue
            x=line.rstrip('\n').split('\t');exact_to_family[x[1].split()[0]]=x[0].split()[0]
    missing=exact_reps-set(exact_to_family)
    if missing:raise SystemExit(f'FAIL family clustering missing exact reps n={len(missing)}')
    fam_candidates=defaultdict(list);fam_exact=defaultdict(set)
    for q in base:
        er=exact[q];fr=exact_to_family[er];fam_candidates[fr].append(q);fam_exact[fr].add(er)

    with open(out/'F3B2_CANDIDATE_FAMILY_MAP.tsv','w') as f:
        f.write('candidate_id\toriginal_query_name\tfamily_rep_id\texact_rep_id\taa_length\taccession\tformal_model_n\tqualifying_model_n\tbest_full_evalue\tbest_domain_i_evalue\n')
        for q in sorted(base,key=lambda q:int(base[q]['serialized_record_ordinal'])):
            r=base[q];er=exact[q];fr=exact_to_family[er]
            f.write('\t'.join([q,r['original_query_name'],fr,er,r['aa_length'],r['accession'],r['formal_model_n'],r['qualifying_model_n'],r['best_full_evalue'],r['best_domain_i_evalue']])+'\n')

    rows=[]
    for fr,members in fam_candidates.items():
        lengths=[int(base[q]['aa_length']) for q in members]
        acc={base[q]['accession'] for q in members if base[q]['accession']!='UNKNOWN'}
        qnames={base[q]['original_query_name'] for q in members}
        longest=max(members,key=lambda q:(int(base[q]['aa_length']),-int(base[q]['serialized_record_ordinal'])))
        rows.append({'family_rep_id':fr,'candidate_instance_n':len(members),'original_query_name_n':len(qnames),'exact_aa_rep_n':len(fam_exact[fr]),'accession_n':len(acc),'multirun_replication':len(acc)>=2,'min_aa_length':min(lengths),'median_aa_length':statistics.median(lengths),'max_aa_length':max(lengths),'longest_candidate_id':longest,'max_formal_model_n':max(int(base[q]['formal_model_n']) for q in members),'max_qualifying_model_n':max(int(base[q]['qualifying_model_n']) for q in members)})
    rows.sort(key=lambda x:(-x['accession_n'],-x['candidate_instance_n'],-x['max_aa_length'],x['family_rep_id']))
    cols=['family_rep_id','candidate_instance_n','original_query_name_n','exact_aa_rep_n','accession_n','multirun_replication','min_aa_length','median_aa_length','max_aa_length','longest_candidate_id','max_formal_model_n','max_qualifying_model_n']
    with open(out/'F3B2_FAMILY_SUMMARY.tsv','w') as f:
        f.write('\t'.join(cols)+'\n')
        for r in rows:f.write('\t'.join(str(r[c]) for c in cols)+'\n')
    with open(out/'F3B2_FAMILY_TOP100.tsv','w') as f:
        f.write('\t'.join(cols)+'\n')
        for r in rows[:100]:f.write('\t'.join(str(r[c]) for c in cols)+'\n')
    gate={'stage':'Stage05F3B2-R3','gate':'PASS_F3B2_R3_CORRECTED_INSTANCE_FAMILY_COMPRESSION','gate_pass':True,'formal_candidate_instance_n':len(base),'exact_aa_unique_n':len(exact_reps),'sequence_family_n':len(rows),'multirun_family_n':sum(r['multirun_replication'] for r in rows),'max_accession_replication_n':max(r['accession_n'] for r in rows),'max_family_candidate_instance_n':max(r['candidate_instance_n'] for r in rows),'family_min_seq_id':lock['family_min_seq_id'],'family_coverage':lock['family_coverage'],'family_cov_mode':lock['family_cov_mode'],'family_cluster_mode':lock['family_cluster_mode'],'family_sensitivity':lock['family_sensitivity'],'known_rt_annotation_used':False,'novelty_screen_performed':False,'family_clustering_role':'operational sequence-family compression only; not proof of biological family identity','candidate_family_map_sha256':sha256_file(out/'F3B2_CANDIDATE_FAMILY_MAP.tsv'),'family_summary_sha256':sha256_file(out/'F3B2_FAMILY_SUMMARY.tsv')}
    with open(out/'STAGE05F3B2_R3_GATE.json','w') as f:json.dump(gate,f,indent=2);f.write('\n')
    print(json.dumps(gate,indent=2))

if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare(sys.argv[2:])
    elif sys.argv[1]=='finalize':finalize(sys.argv[2:])
    else:raise SystemExit('mode must be prepare or finalize')
