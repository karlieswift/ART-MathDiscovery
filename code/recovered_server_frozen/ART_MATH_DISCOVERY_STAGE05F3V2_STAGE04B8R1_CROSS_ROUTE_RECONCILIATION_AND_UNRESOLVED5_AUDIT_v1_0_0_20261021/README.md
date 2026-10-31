# Stage04B8R1 cross-route reconciliation and unresolved-5 audit

Role: analysis-only semantic correction after the completed Stage04B8 closeout.

This stage performs **no new discovery, no structure prediction, no threshold change, and no retuning**.

It addresses two issues:

1. Reconciles the two distinct 54G DELTA_ONLY routes:
   - protein-priority DELTA32 frozen projection -> Stage04B3 clean-RT structure closure;
   - exhaustive architecture-first DELTA_ONLY strict+recurrent route -> Stage04B7 structure closure.
2. Re-audits the five B6 `ABSENT` unitig cases on their original source contigs. `ABSENT` only meant the candidate DNA had no exact match in the Logan unitig file; it did not itself prove that the source-contig ORF was edge-censored.

The output explicitly prevents interpreting Stage04B8's `0/12` as "54G had no clean-RT structural positives overall".
