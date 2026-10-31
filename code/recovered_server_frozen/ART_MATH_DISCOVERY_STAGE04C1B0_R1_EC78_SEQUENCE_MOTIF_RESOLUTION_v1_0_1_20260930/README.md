# Stage04C1B0-R1 EC78 sequence-motif truth-resolution hotfix

This hotfix supersedes the failed Stage04C1B0 v1.0.0 truth-resolution attempt only.
It does not read, recompute, mutate, or rank any Stage04C1 score.

Failure cause in v1.0.0: the resolver required GenBank product text to label an Ec78 RT and nearby ATPase/HNH genes. CP026202.1 annotations do not provide a unique match under that text heuristic, so the run correctly failed closed with 0 candidates.

R1 uses published sequence-level properties of the native Ec78/Eco7 RT: 311 aa and an exact YADD catalytic core. It reuses the already-downloaded GenBank flatfiles from the failed v1.0.0 directory and performs no network fetch.

The resolver fails closed unless exactly one CP026202.1 CDS is 311 aa and contains exact `YADD`. PtuA/PtuB coordinates are not guessed and are not required for the primary cross-control RT-anchor endpoint.
