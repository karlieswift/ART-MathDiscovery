# Stage04B-R1D1 official-protein audit v1.0.1

This is a schema/linkage audit revision, NOT a new mathematical screen and NOT
an automatic override of the v1.0.0 HOLD. Use a new directory; keep old results.
No network calls, BLAST, genome prediction, RT-domain scan, or UG27 changes.

## Run

```
conda activate bioh200
./run_selftest.sh
./run_official_protein_audit.sh /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
./collect_results.sh
```

The audit runner returns 0 only for a complete official-protein linkage audit,
with the previously declared NC_138483.1 genome-only gap. It returns 2 on HOLD.
Run collect_results.sh separately after either outcome; do not delete the HOLD.

Root overrides are supported by the audit script:
--stage04a-root, --r1a-root, --r1b-root, --r1c-root.

## What changed

- Traverse genes[].cds[] and all schema-declared maturePeptide child nodes.
- Keep feature_type and parent_json_path: processed products are not counted as
  independent genes, candidates, or independent discoveries.
- Require explicit identifier links and exact genome/CDS reconstruction for each
  mapped record. Protein sequence always comes from frozen protein.faa.
- Inventory both unlinked CDS FASTA and unlinked protein FASTA records. Missing
  links cannot be bypassed because the number of issues is small.
- Record the exact JSON paths of all sequence identifiers. Return the original
  frozen ZIP for every unresolved accession, plus three fixed diagnostic ZIPs.
- Separate CDS FASTA payload orientation from provider-declared coding direction.
  The historical cds_strand field was not a reliable protein-coding direction.
  New cds_strand is '?' when coding direction is not explicitly declared;
  cds_fasta_payload_strand retains the exact export/reconstruction orientation.
- Compare old vs official AA sequences by an unambiguous coordinate locus; record
  coding-strand disagreements separately rather than conflating them with AA
  identity. Difference patterns are descriptive, not causal explanations.

## Unchanged

All 383 candidate rows, their score/rank/coordinates, +/-10000 bp clipped linear
neighborhoods, canonical source hashes, and the explicit genome-only gap. No
primary retranslation, stop-to-X conversion, candidate dropping, or RT-negative
calls. Circular wraparound, coding-strand inference, genetic-code selection,
and additional ORF prediction are not introduced by this version.

## Outputs

GATE_STAGE04B_R1D1.json
OFFICIAL_NEIGHBORHOOD_PROTEINS.faa
OFFICIAL_NEIGHBORHOOD_PROTEIN_REGISTRY.tsv
CANDIDATE_PROTEIN_LINKS.tsv
ALL_383_CANDIDATES.tsv
R1D0_VS_OFFICIAL_PROTEIN_AUDIT.tsv
REPORT_NODE_TRAVERSAL.tsv
ALL_JSON_FASTA_REFERENCES.tsv
UNLINKED_FASTA_INVENTORY.tsv
ANNOTATION_MAPPING_ISSUES.tsv
RT_RELATED_TEXT_DETAIL.tsv
source_evidence/*.zip and original R1D0 FASTA/registry

Existing source hash locks are checked on every run. Official records whose
protein coding direction remains undeclared can still be sequence-domain input:
that does not authorize any strand-dependent biological claim. Any unresolved
CDS/protein link, unsupported feature, or unresolved old locus keeps the full
linkage gate on HOLD. The genome-only candidate requires a separate ORF route.

## Validation scope

16 local synthetic tests passed. Three use a full 174-genome/383-candidate
synthetic fixture: ordinary input, nested processed products, and a deliberately
unlinked source record. Real results for v1.0.1 have NOT been run here, because
this chat's returned archive does not contain the frozen source ZIPs/R1D0 files.
Do not cite a synthetic self-test as a real biological or source-linkage PASS.

Read docs/RETURNED_V1_0_0_AUDIT.md for the actual returned-data audit.
