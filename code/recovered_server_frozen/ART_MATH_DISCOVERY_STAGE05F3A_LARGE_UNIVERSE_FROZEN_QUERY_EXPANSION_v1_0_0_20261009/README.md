# Stage05F3 LARGE-UNIVERSE — F3A frozen-query expansion

This package formally opens the F3 LARGE-UNIVERSE program.

## F3A0: query-blind source freeze
A 5–10 Gb external raw-DNA universe is built from public Logan metagenomic/virome contigs. Source selection uses SRA metadata only and is completed/frozen before any F2A query sequence is loaded. Selection is balanced by eight predeclared environmental strata using stable accession hashes; individual accessions contributing >250 Mb of eligible contigs are skipped to reduce domination by single runs. Contigs shorter than 300 nt are excluded.

## F3A1: 30 frozen-query expansion
All 30 Stage05F2A novelty-first strict-rescue proteins are searched against the frozen external nucleotide universe using MMseqs translated protein→nucleotide search. Returned target ORFs are reconstructed annotation-independently as maximal stop-to-stop peptides and rescored with the frozen ~45-profile RVT-All HMM collection. F2E homology/longer/HMM rules are unchanged.

Primary endpoint: number of the 26 frozen F2E query-components with an external non-exact longer homolog that retains the frozen formal RVT-All criterion.

## F3B is already method-frozen
`docs/F3B_FROZEN_PROTOCOL.md` freezes the de novo discovery design now. F3B executable code is generated only after F3A0 produces the external FASTA SHA256, and it must reuse exactly that source universe. This prevents source changes after observing F3A outcomes.

Large external FASTA / downloaded Logan contigs are not included in the result collector.
