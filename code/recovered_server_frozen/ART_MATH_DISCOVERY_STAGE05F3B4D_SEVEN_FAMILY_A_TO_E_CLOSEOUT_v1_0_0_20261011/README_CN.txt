Stage05F3B4D — Seven-family A–E computational closeout
Version: v1.0.0 / 2026-10-11

Frozen target:
  all and only the 7 Stage05F3B4C families with:
    strict_published_ART_geometry + recurrent nearest-downstream partner cluster
  No manual family addition/removal after downstream results.

A. Near-context de-replication
  - core architecture = leftmost strict repeat start -> recurrent nearest partner ORF stop
  - exact SHA256 grouping
  - near-context grouping: 31-mer Jaccard >=0.95 and length ratio >=0.95
  - this is de-replication, NOT proof of independent evolutionary lineages

B. Repeat-core conservation
  - all strict arrays from the 7 frozen families
  - no detector retuning
  - one-edit motif-core neighborhoods (Levenshtein <=1)
  - copy/spacer/distance conservation
  - detector motif core != full biological repeat-unit consensus

C. Partner conservation
  - only recurrent strict-array dominant nearest-partner cluster
  - exact AA recurrence, length, orientation, RT-partner gap
  - descriptive Pfam/PHROG/VOG annotation if local resources exist
  - nearest ORF remains a partner candidate, not confirmed dedicated partner

D. Frozen R5C structure extension
  - one architecture-linked representative from each of the 6 target families not already structure-tested
  - selection: maximize min(upstream,downstream), then total context, then candidate_id lexical
  - exact frozen Boltz/R5C settings: Boltz2, server MSA, seed 20261007, recycling=3,
    sampling=200, diffusion=1, no_kernels, serial-exclusive idle GPU; no truncation/fallback
  - STRONG: TM>=0.50 & qcov>=0.50 & tcov>=0.50
  - SUGGESTIVE: TM>=0.40 & qcov>=0.30 & tcov>=0.30

E. Known-system / novelty crosswalk
  - same frozen local known-reference panel construction and MMseqs settings as F3B3B
  - optional user-supplied global FASTA via F3B4D_GLOBAL_REFERENCE_FASTAS
  - optional complete ART95 FASTA via F3B4D_FULL_ART95_FASTA
  - WITHOUT comprehensive global sequence reference + full ART95 sequence catalogue,
    Anthropic-unreported novelty is NOT authorized.

Recommended execution (parallel after prepare):
  1) conda activate art_novelty
     ./run_stage05f3b4d_prepare.sh
  2a) CPU terminal, art_novelty:
     ./run_stage05f3b4d_ace.sh 32
  2b) idle-H200 terminal, art_structure:
     ./run_stage05f3b4d_d_boltz.sh auto
  3) after Boltz, art_novelty:
     ./run_stage05f3b4d_d_foldseek.sh 32
  4) art_novelty:
     ./finalize_stage05f3b4d.sh
     ./collect_stage05f3b4d.sh

Optional deeper partner VOG annotation:
  export F3B4D_RUN_VOG=1

Optional comprehensive references, if/when available:
  export F3B4D_GLOBAL_REFERENCE_FASTAS=/path/a.faa:/path/b.faa
  export F3B4D_FULL_ART95_FASTA=/path/ART95_complete.faa

Claim boundary:
  A family with strict array + recurrent nearest partner + R5C STRONG/SUGGESTIVE is a
  high-priority computational ART-like candidate system. It is not a confirmed ART,
  not a biochemically confirmed RT system, and not demonstrably Anthropic-unreported
  until complete novelty crosswalk and functional validation are finished.
