# Stage05F1 — Annotation-independent RT rescue

Purpose: scan the frozen 3,674 dual-arm candidate neighborhoods directly from nucleotide sequence, without requiring a predicted start codon or an existing gene model, and compare rescued six-frame peptides back to the existing Prodigal/HMM evidence only after rescue.

Primary neighborhood: each frozen dual-arm interval ±10 kb, clipped to the genome. Overlapping neighborhoods are merged per genome for translation efficiency while mapping peptides back to every overlapping candidate.

Outputs include:
- frozen source resolution report
- six-frame peptide registry and FASTA
- candidate↔peptide map
- post-rescue comparison with the existing predicted-ORF registry when available
- RVT-All HMM screen using the existing local panel
- YxDD-like motif support
- rescue prioritization table
- optional MMseqs cross-genome family clustering

No external database download is needed. This can run in parallel with E3H Pfam/PHROG/VOG downloads.
