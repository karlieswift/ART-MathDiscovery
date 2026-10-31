# Stage04B-R1D1: official-protein provenance audit

Version 1.0.0 / 2026-10-05. Python >=3.10, standard library only for the audit;
use the existing bioh200 environment. No downloads or installations occur.

## Why this audit is required

R1D0 generated and hashed 9,560 translated records linked to 382 candidates.
That is an output-generation result, not proof that every translated amino-acid
sequence is the official protein. The previous translation used the default
code, did not validate complete-CDS starts and replaced internal '*' with 'X'.
R1C also guessed strands from simple headers, skipped unsupported headers, and
did not cap right neighborhood endpoints at the genome length.

This package preserves those files, audits them and creates a separate official
protein universe. It does not force the corrected count to remain 9,560.
ALL 383 mathematical candidates and their original scores/ranks are retained.

## Sources and matching

1. Verify the canonical 383-region registry and 240-sequence Stage04A manifest.
2. Verify the R1A manifest and EACH cached ZIP, not just loose unpacked files.
3. Recheck every frozen FASTA file hash and every Datasets genome sequence.
4. Join annotation_report.jsonl genes[].cds[] nucleotide/protein identifiers to
   cds.fna/protein.faa by exact seqId or unique accession-version link.
5. Convert the provider's 1-based ranges to 0-based half-open parts. Verify that
   the oriented/ordered genomic parts reproduce the stored CDS EXACTLY. Explicit
   slippage, mixed strands or ambiguous mapping is recorded as an issue.
6. Rebuild fixed +/-10 kb neighborhoods, clipped to the native record endpoints.
   Keep compound parts rather than assuming one header always gives one interval.
7. Copy the authoritative amino-acid sequence from protein.faa without retranslation.
8. Compare ALL old R1D0 records to the official locus-linked records and report
   identical, different, changed-strand, ambiguous and unresolved cases separately.

No protein is linked by FASTA row order, nearest-name matching, a global best
sequence hit or an RT keyword. No unsupported CDS is silently skipped. A gap in
annotation is not a negative biological call. Unknown 'X' residues supplied by
the official protein are retained and counted; new 'X' residues are not invented.

## Run

```bash
conda activate bioh200
./run_selftest.sh
./run_official_protein_audit.sh /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cat outputs/GATE_STAGE04B_R1D1.json
./collect_results.sh
```

The script uses the exact existing Stage04A, R1A, R1B and R1C root names already
established in this project. Optional --stage04a-root, --r1a-root, --r1b-root and
--r1c-root arguments can specify relocations; contents are still hash-checked.
A rerun archives this package's previous outputs under prior_runs/ rather than
deleting upstream history. No Stage04A/R1/R1C/R1D0 source output is overwritten.

## Important outputs

- ALL_383_CANDIDATES.tsv: every candidate, unchanged mathematical columns.
- OFFICIAL_NEIGHBORHOOD_PROTEINS.faa and OFFICIAL_NEIGHBORHOOD_PROTEIN_REGISTRY.tsv.
- R1D0_VS_OFFICIAL_PROTEIN_AUDIT.tsv: all 9,560 former records reconciled.
- CANDIDATE_PROTEIN_LINKS.tsv and R1C_NEIGHBORHOOD_AUDIT.tsv.
- ANNOTATION_MAPPING_ISSUES.tsv: no silent parse failures.
- RT_RELATED_TEXT_DETAIL.tsv: exact product, protein ID, locus distance.
- PREEXISTING_CDD_ANNOTATIONS.tsv: provider-supplied domain annotation only;
  this is NOT a newly run independent domain search.
- genome_only_rescue/: NC_138483.1 full frozen genome and a fixed-neighborhood
  handoff for later ORF prediction. No gene predictor is run here.

Return the generated RESULTS tar.gz and .sha256. It contains derived sequences
and audit tables, but not the 28 MB cache or the 517 MB hCom2 archive.

## Gate meaning

PASS_OFFICIAL_PROTEIN_AUDIT_WITH_DECLARED_GENOME_ONLY_GAP requires validated
source linkage, all 383 retained candidates, all 174 exact genomes, no silent
CDS mapping issues and no unresolved/ambiguous old-record reconciliation.
Amino-acid differences from the old automatic translation are reported; they
do not invalidate authoritative source retrieval.

A PASS makes the OFFICIAL protein universe ready for an RT-domain scan. It does
not complete that scan, establish ART identity, establish novelty or turn the
unannotated genome into a negative. The genome-only ORF branch remains pending.
Any HOLD includes explicit issue tables; do not raise similarity thresholds,
remove inconvenient candidates or overwrite the previous freeze to get a PASS.

Local tests include a full SYNTHETIC 174-genome / 383-candidate / 9,560-record
structural regression. They are not the actual user's protein comparison.
The actual number of differing translations is unknown until this audit runs.
