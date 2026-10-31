#!/usr/bin/env python3
import csv
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

FORMAL_E = 1e-4
EXPLORE_E = 1e-2


def sha256_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''):
            h.update(b)
    return h.hexdigest()


def main():
    if len(sys.argv) != 9:
        raise SystemExit('usage: collect_instance_hmm.py TBL DOM MAP FROZEN_QNAMES EXTRACT_AUDIT OUTDIR HMM_SHA PROTOCOL_SHA')
    tbl, dom, map_tsv, frozen_qnames, extract_audit, outdir, hmm_sha, protocol_sha = sys.argv[1:9]
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)

    meta={}
    with open(map_tsv) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            q=r['instance_id']
            if q in meta: raise SystemExit(f'FAIL duplicate instance id in map {q}')
            meta[q]=r
    if not meta: raise SystemExit('FAIL empty instance map')

    old=[x.strip() for x in open(frozen_qnames) if x.strip()]
    if len(old)!=10464 or len(set(old))!=10464:
        raise SystemExit('FAIL frozen formal query-name set')
    oldset=set(old)

    full={}
    with open(tbl) as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            x=line.split()
            profile=x[0]; query=x[2]; e=float(x[4])
            if query not in meta: raise SystemExit(f'FAIL unknown query in tblout {query}')
            k=(query,profile)
            if k not in full or e<full[k]: full[k]=e

    ie={}
    with open(dom) as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            x=line.split()
            profile=x[0]; query=x[3]; v=float(x[12])
            if query not in meta: raise SystemExit(f'FAIL unknown query in domtblout {query}')
            k=(query,profile)
            if k not in ie or v<ie[k]: ie[k]=v

    rows=[]; qrows=defaultdict(list)
    for k,e in full.items():
        if k not in ie: continue
        d=ie[k]
        if e<=FORMAL_E and d<=FORMAL_E: cls='FORMAL'
        elif e<=EXPLORE_E and d<=EXPLORE_E: cls='EXPLORATORY'
        else: continue
        q,p=k
        rows.append((q,meta[q]['original_query_name'],p,e,d,cls))
        qrows[q].append((p,e,d,cls))

    rows.sort(key=lambda z:(int(meta[z[0]]['serialized_record_ordinal']),z[2]))
    with gzip.open(out/'STAGE05F3B1R3_INSTANCE_RVT_PAIR_HITS.tsv.gz','wt') as f:
        f.write('query_id\toriginal_query_name\tprofile\tfull_evalue\tbest_domain_i_evalue\tclass\n')
        for r in rows:
            f.write('%s\t%s\t%s\t%.8g\t%.8g\t%s\n'%r)

    formal=[]; exploratory=[]
    with open(out/'STAGE05F3B1R3_INSTANCE_CANDIDATE_SUMMARY.tsv','w') as f:
        f.write('query_id\toriginal_query_name\tclass\tmodel_n\tbest_full_evalue\tbest_domain_i_evalue\n')
        ordered=sorted(qrows,key=lambda q:int(meta[q]['serialized_record_ordinal']))
        for q in ordered:
            rr=qrows[q]
            cls='FORMAL' if any(x[3]=='FORMAL' for x in rr) else 'EXPLORATORY'
            (formal if cls=='FORMAL' else exploratory).append(q)
            f.write(f"{q}\t{meta[q]['original_query_name']}\t{cls}\t{len(set(x[0] for x in rr))}\t{min(x[1] for x in rr):.8g}\t{min(x[2] for x in rr):.8g}\n")

    (out/'STAGE05F3B1R3_FORMAL_INSTANCE_IDS.txt').write_text(''.join(q+'\n' for q in formal))
    (out/'STAGE05F3B1R3_EXPLORATORY_ONLY_INSTANCE_IDS.txt').write_text(''.join(q+'\n' for q in exploratory))

    recovered_qnames={meta[q]['original_query_name'] for q in formal}
    missing=oldset-recovered_qnames
    extra=recovered_qnames-oldset
    if missing or extra:
        raise SystemExit(f'FAIL corrected formal-instance rerun does not recover frozen formal query-name set: missing={len(missing)} extra={len(extra)} missing_examples={sorted(missing)[:20]}')

    formal_by_name=defaultdict(int)
    for q in formal: formal_by_name[meta[q]['original_query_name']]+=1
    a=json.load(open(extract_audit))

    gate={
      'stage':'Stage05F3B1R3',
      'gate':'PASS_F3B1R3_QUERY_NAMESPACE_RESCUE_COMPLETE',
      'gate_pass':True,
      'old_f3b1_formal_query_name_n':10464,
      'corrected_formal_instance_n':len(formal),
      'corrected_formal_query_name_n':len(recovered_qnames),
      'exploratory_only_instance_n':len(exploratory),
      'qualifying_pair_n':len(rows),
      'formal_query_name_with_multiple_formal_instances_n':sum(v>1 for v in formal_by_name.values()),
      'max_formal_instance_n_per_original_query_name':max(formal_by_name.values()),
      'selected_namespace_instance_n':a['selected_namespace_instance_n'],
      'formal_query_name_with_duplicate_input_instances_n':a['formal_query_name_with_duplicate_instances_n'],
      'formal_full_evalue_threshold':FORMAL_E,
      'formal_domain_i_evalue_threshold':FORMAL_E,
      'exploratory_full_evalue_threshold':EXPLORE_E,
      'exploratory_domain_i_evalue_threshold':EXPLORE_E,
      'rvt_hmm_sha256':hmm_sha,
      'protocol_sha256':protocol_sha,
      'annotation_used_for_candidate_generation':False,
      'scientific_threshold_changed':False,
      'old_10464_interpretation':'unique HMM query-name count after namespace collapse; not a guaranteed peptide-instance count',
      'corrected_candidate_identity':'unique peptide instances identified by global HMM-input FASTA record ordinal and validated against original MMseqs ORF DB',
      'tblout_sha256':sha256_file(tbl),
      'domtblout_sha256':sha256_file(dom),
      'claim_boundary':'Query-namespace rescue only; the frozen F3B1 RVT-profile screen and thresholds are unchanged.'
    }
    with open(out/'STAGE05F3B1R3_GATE.json','w') as f:
        json.dump(gate,f,indent=2);f.write('\n')
    print(json.dumps(gate,indent=2))

if __name__=='__main__': main()
