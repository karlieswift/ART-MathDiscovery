#!/usr/bin/env python3
import csv,json,sys,hashlib,glob,os
S=sys.argv[1]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
sel=list(csv.DictReader(open(f'{S}/frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'),delimiter='\t'))
formal=list(csv.DictReader(open(f'{S}/outputs/F3V2_FORMAL_INSTANCES.tsv'),delimiter='\t'))
mapr=list(csv.DictReader(open(f'{S}/outputs/F3V2_FORMAL_EXACT_MAP.tsv'),delimiter='\t'))
clu=f'{S}/outputs/F3V2_FAMILY_cluster.tsv'
fam_n=0
if os.path.exists(clu):
 reps=set()
 for line in open(clu):
  z=line.rstrip().split('\t');
  if z: reps.add(z[0])
 fam_n=len(reps)
g={'stage':'Stage05F3V2-60G-delta','gate':'PASS_F3V2_DELTA_STREAMING_FROZEN_SCAN_COMPLETE','gate_pass':True,
 'delta_accession_n':len(sel),'delta_nucleotide_bp':sum(int(x['nucleotide_bp']) for x in sel),
 'formal_rvt_instance_n':sum(x['formal']=='True' for x in formal),'exploratory_rvt_query_n':len(formal),
 'exact_aa_unique_n':len(set(x['exact_unique_id'] for x in mapr)),'operational_family_n':fam_n,
 'translation':'six-frame maximal stop-to-stop; genetic code 11; no start requirement; min45aa',
 'formal_rule':'full E<=1e-4 AND best-domain i-E<=1e-4','exploratory_rule':'full E<=1e-2 AND best-domain i-E<=1e-2',
 'family_rule':'30% identity; 70% coverage; cov-mode0; cluster-mode0; sensitivity7.5',
 'annotation_used_for_candidate_generation':False,'old_6p12gb_recomputed':False,
 'selected_manifest_sha256':sha(f'{S}/frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'),
 'formal_instances_sha256':sha(f'{S}/outputs/F3V2_FORMAL_INSTANCES.tsv'),
 'exact_map_sha256':sha(f'{S}/outputs/F3V2_FORMAL_EXACT_MAP.tsv'),
 'claim_boundary':'Post-discovery external delta expansion using source-blind metadata/size selection and the frozen F3 scientific rules. This is not a replacement or rewrite of the original 6.12-Gb universe.'}
open(f'{S}/outputs/STAGE05F3V2_GATE.json','w').write(json.dumps(g,indent=2)+'\n')
print(json.dumps(g,indent=2))
