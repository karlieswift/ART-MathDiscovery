# Stage05C2B: unambiguous ten-locus ledger, cassette variation, and RNA hypotheses

## Scientific purpose
Continue towards credible novel ART-like systems without relabelling known-family
relatives, mathematical misses, or sequence folding as biological discovery.
This is a NEW post-score focused follow-up, not a rerun of 919 genomes.

The two whole-genome arms were parallel:
- Frozen DNA mathematics: 919 genomes -> 3,674 candidate regions -> fixed +/-10 kb neighborhoods.
- Uniform Prodigal + RT profiles: all 919 genomes -> 10 RT-core GA loci.
Joining the arms gives 7 covered loci and 3 not covered. Within the 7 are 5 close
three-component ART-reference relatives, 1 distant co-occurrence, and 1 RT-only case.
OR770614.1 is one of the outside 3, not an eleventh locus. It is a close reference
relative, so all six close reference-related cassettes consist of captured 5 plus
this one baseline-only comparator. None is declared a new ART family here.

IMPORTANT: C1 exact-repeat+RT count 5 is NOT the same set as C2A close-cassette count
5. The former includes NC_016163.1, the latter substitutes OR455461.1 (the C1
approximate-repeat case). See precomputed/TEN_LOCUS_UNAMBIGUOUS_LEDGER.tsv.

## Already performed (bundled, not fictional future output)
- Ten-locus orthogonal count/identity ledger based on C1 and C2A source records.
- 45 global component alignments: 15 pairs among the SIX scoped close cassettes,
  separately for RT, partner-like ORF, and COMPLETE RT-oriented immediate upstream
  DNA. These are new continuous descriptive comparisons, not a family threshold,
  not a global database search, and do not overwrite prior local alignments.
- Explicit substitution and indel-column evidence for every upstream pair.
- Frozen 2,396 RNA coordinate/orientation occurrences, 2,324 exact sequence targets:
  all ten RT-centered contexts and three registered reference cores. Entire
  intergenic windows are relative to the frozen Prodigal calls, not transcript
  truth. Both directions retained. Closest-upstream length 1 nt for OK490448.1
  does NOT remove its other RT-context intergenic sequences.

## Still to execute on the server
1. Fetch the official Anthropic technical report. Save bytes/hash, extract URI links,
   embedded files and exact accession mention locations. A mention is not confirmed
   catalogue membership; a missing mention never proves novelty. The full author
   family catalogue has NOT been acquired or compared. Do not use a third-party
   reproduction repository as an author catalogue.
2. RNAfold 2.7.2 MFE calculations for all frozen targets, each with 99 independent
   seeded randomized Euler-trail surrogates. Dinucleotide counts and endpoints are
   exact. This sampler is not asserted uniform; the reported lower-tail proportion
   is NOT a calibrated p-value, not FDR and not an ART-positive criterion.

These are computational hypotheses, not wet-lab constructs. No expression,
RT activity, partner function, structural conservation or new family is proved.

## Setup (new environment; existing environments and UG27 unchanged)

    conda create -n art_rna_review -y --override-channels --strict-channel-priority \
      -c conda-forge -c bioconda \
      python=3.11 biopython=1.86 viennarna=2.7.2 pymupdf
    conda activate art_rna_review
    ./run_selftest.sh
    ./run_stage05c2b.sh 4

All analysis input DNA/protein is bundled; no NCBI download or HMMER/Prodigal rerun.
An unavailable author report leaves an explicit source HOLD; it does not prevent
RNA characterization. If a verified local author PDF is already present:

    ./fetch_author_sources.sh /absolute/path/technical_report.pdf

To replay the new precomputed comparisons independently of RNAfold:

    ./reproduce_precomputed.sh

Restart the same run command to reuse validated completed target caches. A changed
binary, Python runtime, policy or target input refuses to mix with the old cache.
No whole-genome scoring, original candidate edits, parameter retuning, or UG27 writes.

## Results and return package

    cat outputs/GATE_STAGE05C2B.json
    cat outputs/author_source/AUTHOR_SOURCE_GATE.json
    cat outputs/RNA_FOLD_GATE.json
    ./collect_results.sh

Upload ART_STAGE05C2B_TEN_LOCUS_RESULTS_*.tar.gz and .sha256. Raw RNAfold input,
output and cache metadata are retained. The official PDF, if acquired, is included
for actual source review; full extracted article text is excluded. A technically
complete RNA execution does not change the outstanding full-catalogue novelty HOLD.

## Scope and provenance
Original C1 archive sha256:
e0f2e8e5049d2e7ec9526a0c1d1a8cae08d3e32bfa2866a66e5680f8711772fa
Original C2A artifact sha256:
bc0addd8ee17093f4a61895843c25179a0df70ad88a14335f93f0d591f48fdb2

All 3,674 math candidates remain upstream. The ten-locus follow-up is scoped by
profile evidence; candidates without these profiles are NOT biologically negative.
The six close-cassette subset is post-hoc and cannot estimate accuracy or novelty.
