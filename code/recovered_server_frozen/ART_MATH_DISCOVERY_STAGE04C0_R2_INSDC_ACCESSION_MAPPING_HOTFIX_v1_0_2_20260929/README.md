# Stage04C0-R2 INSDC accession-mapping hotfix v1.0.2

This package supersedes the failed Stage04C0-R1 sequence acquisition attempt.

Root cause: DDBJ GetEntry serves INSDC/GenBank/DDBJ/EMBL accessions, whereas two C0 rows used NCBI RefSeq `NC_` accessions. `NC_016856.1` and `NC_012731.1` therefore returned tiny non-FASTA responses from DDBJ. The sequences themselves are valid RefSeq records; the provider/accession namespace was the problem.

Frozen source mapping:
- EC86_ECO1: `CP001509.3 -> CP001509.3`
- SEN2_14028S: `NC_016856.1 -> CP001363.1`
- EC78_ECO7: `CP026202.1 -> CP026202.1`
- DRT2_KPN: `NC_012731.1 -> AP006725.1`

The two mapped source records are the corresponding complete GenBank/DDBJ records for the same genomes. No ACNV or recurrence scoring is performed here.

Successful final gate: `PASS_STAGE04C0_R2_INSDC_ACCESSION_MAPPING_FIXED`.

C1 should use the R2 root, not the failed R1 root.
