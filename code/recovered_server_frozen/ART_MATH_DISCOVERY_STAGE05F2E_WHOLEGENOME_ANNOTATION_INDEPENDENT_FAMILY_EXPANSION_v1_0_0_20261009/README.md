# Stage05F2E v1.0.0

Whole-genome annotation-independent family expansion for the 30 frozen novelty-first
Stage05F2A strict rescue sequences.

Why this stage exists
---------------------
Stage05F2D searched only the Stage05F1 +/-10 kb candidate-neighborhood six-frame
universe. That is the correct discovery-context universe, but it cannot answer
whether a longer homolog of a rescued family occurs elsewhere in the same frozen
919 genomes. F2E removes that remaining representation limitation without changing
the frozen DNA detector, candidate set, or query set.

F2E translates all 919 frozen source genomes in six frames as maximal stop-to-stop
peptides (no start-codon requirement, minimum 45 aa), searches all 30 frozen F2A
queries against that whole-genome peptide universe, and runs the already frozen
RVT-All profile collection only on returned homolog sequences. Existing Prodigal
ORFs and the frozen 3,674 candidate +/-10 kb neighborhoods are consulted post hoc.

Required inputs
---------------
1. Stage05F2A live directory or collected result tar.
2. Stage05F2D live directory or collected result tar.
3. ART_MathDiscovery project root containing the frozen Stage05A/05B source files,
   Stage05C1 ORF registry if available, and frozen RVT-All HMM.
4. `mmseqs` and `hmmscan` in PATH.

Run
---
    ./run_stage05f2e.sh F2A_PATH F2D_PATH PROJECT_ROOT [THREADS]

The output gate is `outputs/STAGE05F2E_GATE.json`.
