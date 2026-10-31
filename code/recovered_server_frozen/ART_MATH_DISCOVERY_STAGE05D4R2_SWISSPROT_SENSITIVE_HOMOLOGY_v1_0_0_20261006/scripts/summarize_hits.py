#!/usr/bin/env python3
from pathlib import Path
import csv,re,json,sys,collections
qfa,tfa,ctab,hits,out=map(Path,sys.argv[1:6]); out.mkdir(parents=True,exist_ok=True)
def headers(p):
    d={}
    with p.open(errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                s=line[1:].rstrip('\n'); key=s.split()[0]; d[key]=s
    return d
th=headers(tfa)
# query headers are AA_<sha>; cluster representatives table maps bare sha
rtpat=re.compile(r'reverse transcript|RNA[- ]directed DNA polymerase|retron|diversity[- ]generating retro|group II intron|maturase|CRISPR.*reverse transcript',re.I)
known_bad=re.compile(r'virion structural|capsid|tail protein|terminase|portal protein',re.I)
rows=[]
if hits.exists() and hits.stat().st_size:
    with hits.open() as f:
        rd=csv.reader(f,delimiter='\t')
        for r in rd:
            if len(r)<12: continue
            q,t=r[0],r[1]; desc=th.get(t,'')
            rows.append({'query':q,'target':t,'fident':float(r[2]),'alnlen':int(r[3]),'qstart':int(r[4]),'qend':int(r[5]),'tstart':int(r[6]),'tend':int(r[7]),'qcov':float(r[8]),'tcov':float(r[9]),'evalue':float(r[10]),'bits':float(r[11]),'target_description':desc,'rt_text':bool(rtpat.search(desc)),'structural_text':bool(known_bad.search(desc))})
# all query ids
qs=[]
with qfa.open() as f:
    for line in f:
        if line.startswith('>'): qs.append(line[1:].split()[0])
byq=collections.defaultdict(list)
for r in rows: byq[r['query']].append(r)
summary=[]
for q in qs:
    arr=sorted(byq.get(q,[]), key=lambda x:(x['evalue'],-x['bits']))
    best=arr[0] if arr else None
    rt=next((x for x in arr if x['rt_text']),None)
    summary.append({'query':q,'n_hits':len(arr),'best_target':best['target'] if best else '','best_evalue':best['evalue'] if best else '','best_bits':best['bits'] if best else '','best_qcov':best['qcov'] if best else '','best_tcov':best['tcov'] if best else '','best_description':best['target_description'] if best else '','any_rt_text_hit':rt is not None,'best_rt_target':rt['target'] if rt else '','best_rt_evalue':rt['evalue'] if rt else '','best_rt_qcov':rt['qcov'] if rt else '','best_rt_tcov':rt['tcov'] if rt else '','best_rt_description':rt['target_description'] if rt else '','best_is_structural_text':bool(best and best['structural_text'])})
with (out/'STAGE05D4R2_QUERY_SUMMARY.tsv').open('w',newline='') as f:
    fields=list(summary[0].keys()); w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(summary)
# candidate/family mapping using NR90 cluster table representative -> members/candidate keys
repmap=collections.defaultdict(lambda:{'candidate_keys':set(),'strata':set(),'members':0})
with ctab.open() as f:
    rd=csv.DictReader(f,delimiter='\t')
    for r in rd:
        rep='AA_'+r['representative_aa_sha256']; x=repmap[rep]; x['members']+=1
        x['candidate_keys'].update(filter(None,r['member_candidate_keys'].split(';')))
        x['strata'].update(filter(None,r['member_strata'].split(';')))
tri=[]
for s in summary:
    x=repmap.get(s['query'],{'candidate_keys':set(),'strata':set(),'members':0})
    # triage label is descriptive, not filtering
    if s['any_rt_text_hit']: label='CURATED_RT_TEXT_HOMOLOGY_TRIAGE'
    elif s['n_hits']==0: label='NO_SWISSPROT_HIT_NOT_NEGATIVE'
    elif s['best_is_structural_text']: label='CURATED_STRUCTURAL_PROTEIN_HOMOLOGY_TRIAGE'
    else: label='CURATED_NON_RT_HOMOLOGY_TRIAGE'
    tri.append({**s,'nr90_member_n':x['members'],'candidate_key_n':len(x['candidate_keys']),'candidate_keys':';'.join(sorted(x['candidate_keys'])),'strata':';'.join(sorted(x['strata'])),'triage_label':label})
with (out/'STAGE05D4R2_REMOTE_HOMOLOGY_TRIAGE.tsv').open('w',newline='') as f:
    fields=list(tri[0].keys()); w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(tri)
counts=collections.Counter(x['triage_label'] for x in tri)
report={'n_queries':len(qs),'n_queries_with_any_hit':sum(x['n_hits']>0 for x in tri),'n_queries_with_rt_text_hit':sum(x['any_rt_text_hit'] for x in tri),'n_no_hit':sum(x['n_hits']==0 for x in tri),'triage_counts':dict(counts)}
(out/'STAGE05D4R2_SUMMARY.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
