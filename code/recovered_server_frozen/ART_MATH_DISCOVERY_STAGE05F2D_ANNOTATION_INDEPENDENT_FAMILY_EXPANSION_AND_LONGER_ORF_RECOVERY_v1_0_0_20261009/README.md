# Stage05F2D v1.0.0

Annotation-independent family expansion and longer-ORF recovery triage for all
30 Stage05F2A novelty-first strict rescue sequences.

Required inputs:
1. Stage05F1 live directory or result tar containing the full six-frame outputs.
2. Stage05F2A live directory or result tar containing the 30 novelty-first rows.
3. Optional Stage05F2C-B live directory/result tar; use NONE to omit.
4. MMseqs2 executable in PATH.

The stage searches the 30 frozen queries against the 380,969 unique Stage05F1
six-frame peptides, maps homologous unique sequences back to all 498,879 genomic
instances, and identifies replicated/longer annotation-independent homologs.

No new HMM search or structure prediction is performed. Existing F1 RVT-profile
status is carried forward for homolog targets.
