F3B2 R4 source-accession metadata rescue
========================================

Purpose
-------
R3 successfully rescued candidate identity at the peptide-instance level, but its
source accession parser incorrectly searched SRA/ENA run IDs inside MMseqs ORF
headers. MMseqs extractorfs ORF headers encode the parent input-database key and
ORF coordinates; therefore R3 reported accession_n=0 for all families.

R4 changes metadata/provenance only:
- no HMM rerun
- no candidate membership change
- no exact-AA change
- no family reclustering
- no threshold change
- no novelty/known-RT annotation

R4 follows the parent contig key from each ORF header back to the original
F3B1 genomeDB_h header, validates the accession against the frozen 131-accession
F3A manifest, then recomputes only accession/run replication summaries.
