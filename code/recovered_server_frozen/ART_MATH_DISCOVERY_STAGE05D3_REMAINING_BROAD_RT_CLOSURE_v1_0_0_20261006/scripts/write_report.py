#!/usr/bin/env python3
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs'
a=list(csv.DictReader(open(OUT/'REMAINING_BROAD_RT_HIT_AUDIT.tsv'),delimiter='\t'))
o={r['accession']:r for r in csv.DictReader(open(OUT/'REMAINING_OFFICIAL_ANNOTATION_AUDIT.tsv'),delimiter='\t')}
lines=['# Stage05D3 remaining broad-RT closure','', 'All targets were frozen before Stage05D3; no DNA rescoring/reranking was performed.','']
for r in a:
    x=o.get(r['accession'],{})
    lines += [f"## {r['accession']} / {r['rt_record_key']}",
      f"- Stage05D1 model: {r['myrt_model']}; formal full E={r.get('formal_full_E','')}; domain iE={r.get('formal_domain_iE','')}.",
      f"- ORF length: {r.get('observed_aa_len','')} aa; strand {r['rt_strand']}; frozen candidate distance(s): {r['pre_d3_candidate_distance_bp']} bp; strata {r['pre_d3_strata']}.",
      f"- Motif scan: {r.get('motif_scan','')}.",
      f"- Official annotation: {x.get('status','')}; product={x.get('official_product','')}; protein_id={x.get('official_protein_id','')}; exact translation={x.get('exact_translation_match','')}.",
      '- Interpretation: broad RT-family HMM evidence is not sufficient for RT or ART identity; architecture/annotation conflicts are used only for biological triage.', '']
(OUT/'STAGE05D3_REPORT.md').write_text('\n'.join(lines)+'\n')
