#!/usr/bin/env python3
import argparse,csv,json,re
ap=argparse.ArgumentParser(); ap.add_argument('manifest'); args=ap.parse_args()
rows=list(csv.DictReader(open(args.manifest),delimiter='\t'))
exp={f'ART_{i:02d}' for i in range(1,96)}; ids=[r.get('art_id','').strip() for r in rows if r.get('art_id','').strip()]
seqs=[re.sub('[^A-Z]','',r.get('rt_aa','').upper().replace('*','')) for r in rows if r.get('rt_aa','').strip()]
out={'rows':len(rows),'unique_art_ids':len(set(ids)),'missing':sorted(exp-set(ids)),'duplicates':len(ids)-len(set(ids)),'aa_present':sum(bool(s) for s in seqs),'yxdd':sum(bool(re.search(r'Y.DD',s)) for s in seqs),'ready_for_exact_95_catalogue_claim':len(rows)==95 and set(ids)==exp and len(ids)==95 and sum(bool(s) for s in seqs)==95}
print(json.dumps(out,indent=2,sort_keys=True))
