# Stage05C1 - post-freeze new919 RT / ORF / repeat architecture evidence

This is a NEW biological evidence run, not a Stage05B rerun or a replacement
of the old 240-genome pilot. The goal is to develop credible novel ART-like
system leads for later family/novelty review and wet-lab prioritization.
A protein profile hit, a mathematical rank, and an ART identity are different
claims. All 3,674 frozen mathematical candidates remain in the main output.

## Run on the existing server environment

No internet access, new NCBI download, model update, or new package install is
needed. Activate `art_rt_evidence` (Python >=3.9, Prodigal 2.6.3, HMMER 3.4).
The source and score roots must be the actual existing run directories, not
only the small returned archives. The Stage05B results archive must remain in
its original Stage05B root:
`ART_STAGE05B_NEW919_SCORE_RESULTS_20261006_081746.tar.gz`.
The program verifies its expected SHA256 before biology.

```bash
conda activate art_rt_evidence
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
A05="$P/ART_MATH_DISCOVERY_STAGE05A_NEW_SOURCE_DISCOVERY_FREEZE_v1_0_0_20261005"
B05="$P/ART_MATH_DISCOVERY_STAGE05B_NEW919_FROZEN_TWO_ARM_SCAN_v1_0_0_20261005"
./run_selftest.sh && ./run_stage05c1.sh "$A05" "$B05" 4
```

Use the same command after interruption. Completed genome predictions and
completed candidate repeat analyses are reused ONLY after signature and file
hash verification. Partly finished jobs are recomputed. Do not launch another
writer in this root. Do not modify policy to make a failed job pass.
The program has a write lock and bounded job submission. A failure stops new
submission; already running tasks may finish before shutdown.

## Scope and sequence of computation

1. Check the release files, original result archive, source manifest, all 919
   actual source FASTAs, and ALL four per-genome mathematical score artifacts
   against the existing completion records. This reads, but does not recompute,
   the large window tables. Inspect the same source integrity again at closeout.
2. Run the fixed exact and substitution-tolerant repeat rules on three previously
   registered reference contexts. These are developmental/descriptive capability
   checks, not held-out estimates or gates forcing a favorable biological answer.
3. Run Prodigal meta on all 919 complete genomes, including 91 with no mathematical
   dual-arm candidate. Validate the CDS nucleotide sequence, strand, GFF and
   protein identifiers. Record model-selected translation table and partial flags.
   Remove only a single terminal stop marker; never change internal stops to X.
   Noncanonical or empty proteins produce HOLD, not omission or negative calls.
4. Freeze and scan ALL whole-genome predicted proteins plus six registered RT
   reference proteins, after exact-AA computational deduplication. Restore every
   protein locus and candidate link. Six RT references include three non-array
   RT systems; they are NOT RT-negative controls. The number of RT loci outside
   mathematical neighborhoods is retained, not made into retroactive candidates.
5. Inspect EVERY candidate's fixed +/-10 kb linear, clipped neighborhood for
   exact direct-repeat cores AND a separate ungapped substitution sensitivity.
   No RT result, mathematical score or rank is passed to repeat discovery.
6. Reconstruct the matching k-mer pairs at each frozen best recurrence window/
   k/lag and merge overlapping match positions into paired sequence blocks.
   This is not recalculation of the detector or proof that match count equals
   repeat-copy count. Export the actual sequences and paired coordinates.
7. Combine RT evidence, repeat copies, predicted CDS masks, physical distances,
   strand-relative context, spacer-sequence hashes/exact diversity and nearby ORF inventory into all-candidate evidence
   tables and JSONL cards. Nearby ORFs are not called proven partners.

## Models and limits

The THREE Pfam HMMs are copied byte-for-byte from the prior returned R1D2 panel:
RVT_1/PF00078.33 and RVT_2/PF07727.20 are core evidence; RVT_3/PF13456.13 is
auxiliary RT-like evidence. `--cut_ga` requires model-specific sequence AND
domain gathering scores. A separate declared sensitivity uses full E <=0.001,
independent domain E <=0.001 and model coverage >=0.35. Both scans use seed 42
and fixed `-Z 100000`, as in R1E; this is the declared E-value normalization,
NOT a claim that the actual target count equals 100,000. No-hit is not absence.
No catalytic activity is established by profile membership.

Exact-repeat description reuses the previous exact-core algorithm unchanged:
12-nt seed, >=3 nonoverlapping copies, common-core extension up to 64 nt,
adjacent start distance <=1024 bp. All complexity classes are retained.

The NEW post-score exploratory ungapped sensitivity has fixed core lengths
16/24/32/48, sampled exact 8-mer seeds at 0/middle/end offsets, <=20% substitutions
to an actual sequence template, >=3 nonoverlapping copies, and the same
1024-bp adjacent-start bound. This is NOT an exhaustive approximate aligner:
indels, unsampled-seed differences, inverse repeats and structure-only arrays
can be missed. It is NOT an RNA folding or expression analysis. Its output is
not calibrated to a false-discovery rate. Correlated/contained motif descriptions
may still refer to the same array; counts are NOT independent systems.

Both descriptions compare repeat copies AND the entire group span with predicted
CDS intervals. No predicted overlap means prediction-relative evidence only,
not proven noncoding RNA. Official NCBI annotations are NOT imported in this
stage: these are new de novo predicted proteins, not official protein records.
The fixed linear neighborhoods do not wrap around genome origins; edge clipping
is explicitly flagged. Near homology, known-family exclusion, complete domain/
partner architecture, indel-aware and RNA structural analysis remain later work.

The bounded substitution sensitivity did NOT recover the registered Machias
core in the local reference check either: 2/3 reference cores for exact, 2/3 for
this new sensitivity. It is not presented as a fix that solved every known array.
This limitation is retained; absence is not used to exclude any candidate.

## Main outputs

- `outputs/GATE_STAGE05C1.json`, `SUMMARY_STAGE05C1.json`
- `ALL_3674_SYSTEM_EVIDENCE.tsv`: every original field and row order retained.
- `ALL3674_EVIDENCE_CARDS.jsonl`: every candidate, including unresolved/no-hit.
- `ALL919_RT_PROFILE_LOCI.tsv`: core/auxiliary loci inside AND outside neighborhoods.
- `RT_ALL_REPORTED_DOMAINS.tsv`: raw reported and qualified evidence tiers separated.
- `REPEAT_GROUPS_ALL_CANDIDATES.tsv`: all exact and substitution group descriptions.
- `RT_REPEAT_SYSTEM_CONTEXT.tsv`, `RT_NEIGHBOR_ORF_INVENTORY.tsv`: context, not function.
- `RECURRENCE_MATCH_BLOCK_WITNESSES.tsv`: paired sequence blocks behind k-mer counts.
- `FULL919_GENE_REGISTRY.tsv`, `CANDIDATE_PREDICTED_ORF_LINKS.tsv`.
- `REFERENCE_RT_CALIBRATION.tsv`, `REFERENCE_REPEAT_CALIBRATION.tsv`.

`PASS_ALL3674_POSTFREEZE_RT_ORF_REPEAT_EVIDENCE_WITH_SCOPE_LIMITS` means computation
and linkage are complete; it never means novel ART or wet-lab readiness proven.

## Return results

After completion or after an error has stopped the job:

```bash
./collect_results.sh
```

Return `ART_STAGE05C1_NEW919_BIOLOGY_RESULTS_*.tar.gz` and its `.sha256`.
This includes the entire amino-acid target set, source mapping and evidence
for all candidates. For all qualifying core/auxiliary RT-profile source genomes,
it also copies the original genome and raw Prodigal outputs for direct review.
Other whole-genome FASTAs, bulk original mathematical window tables, all raw
prediction directories and intermediate repeat checkpoints remain on the
server. They are NOT silently described as independently rechecked from the
small returned package. Keep them until the project is closed.

## Local validation identity

46 local unit tests passed. The parser was exercised against actual prior returned files
for 240 genomes / 24,321 predicted records and 16 GA + 30 sensitivity domain
rows. The two-genome integration/continuation/collection fixture uses test stubs
replaying prior real outputs; it is NOT a real HMMER/Prodigal execution.
No new919 biological results are precomputed or bundled as study results.
No UG27, pilot model, threshold, rank, score or account files are touched.
