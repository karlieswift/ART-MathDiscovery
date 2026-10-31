# Stage05F1R1 strict rescue-instance audit v1.0.0

Purpose: deconvolve Stage05F1 AA-level rescue calls at peptide-instance level before any claim that a gene caller missed an RT. No new translation, HMM search, detector scoring, threshold retuning, or structure prediction is performed.

The 41 Stage05F1 `FORMAL_PROFILE_GENE_CALLER_RESCUE` rows are unique-AA classes, not automatically 41 independent RT loci.
