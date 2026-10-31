# Stage05F2E v1.0.2 - exact-recovery hotfix

Whole-genome annotation-independent family expansion for the 30 frozen novelty-first
Stage05F2A strict rescue sequences.

## Why v1.0.2 exists

The v1.0.0 run stopped at:

    exact whole-genome query recovery 5/30

That check was implemented incorrectly: it inferred whether a frozen F2A peptide
exists in the whole-genome six-frame materialization from whether MMseqs returned
an exact hit. MMseqs is a heuristic similarity-search engine; short or low-complexity
phage peptides can be affected by the prefilter and residue-masking path. A returned-hit
list is therefore not a valid generation-integrity invariant.

v1.0.2 fixes the implementation before any F2E biological result is interpreted:

1. Replays every frozen F2A instance coordinate directly from the frozen source FASTA.
2. Audits exact query recovery directly by AA SHA256 in the whole-genome peptide registry.
3. Runs the same frozen family search with MMseqs residue masking disabled (`--mask 0`).
4. If an exact target verified by direct hash is still absent from the MMseqs primary result,
   only that frozen query is rerun with `--exhaustive-search 1`; results are merged deterministically.
5. Records coordinate replay, direct-hash recovery, primary MMseqs recovery, and final MMseqs
   recovery separately in `STAGE05F2E_EXACT_QUERY_RECOVERY_AUDIT.tsv`.

The 30 queries, 919 genomes, translation rule, homology thresholds, RVT-All thresholds,
DNA scores/ranks, and candidate membership are unchanged. This is an implementation
hotfix, not biological retuning.

## Reuse of the failed v1.0.0 run

If the failed v1.0.0 directory still contains complete files:

    work/STAGE05F2E_WHOLE919_SIXFRAME_INSTANCES.faa
    work/STAGE05F2E_WHOLE919_SIXFRAME_INSTANCES.tsv

`run_stage05f2e.sh` automatically hard-links or copies them instead of retranslating all
171,397,507 bp. They are then re-audited by direct AA hashes before reuse.

## Run

    ./run_stage05f2e.sh F2A_PATH F2D_PATH PROJECT_ROOT [THREADS]

The final gate requires 30/30 coordinate replay, 30/30 direct whole-genome AA-hash
recovery, and 30/30 final MMseqs exact-target recovery.

## v1.0.2 translation hotfix

F2E v1.0.0/v1.0.1 accidentally defined only 62/64 codons in its local translation dictionary, omitting AAG->K and GAG->E. Stage05F1 used the complete 64-codon table. This caused coordinate replay mismatches (25/30 frozen queries in the observed v1.0.1 run) and makes the old whole-genome peptide materialization unsuitable for reuse. v1.0.2 restores the exact F1 codon dictionary, locks its SHA256, strengthens selftests, and forces clean regeneration from the frozen 919-genome FASTAs. Biological/query thresholds are unchanged.
