# ART MathDiscovery Stage05D1 — novelty-first mining

Purpose: analyze the **frozen 999 candidate subset** that is most relevant to discovery beyond the current RVT_1/RVT_2/RVT_3 panel:

- 987 `REPEAT_CONTEXT_RT_NOT_DETECTED_BY_THIS_PANEL`
- 12 `RT_LIKE_AUXILIARY_PROFILE_CONTEXT`

This stage does **not** alter Stage05B scores/ranks or Stage05C1 labels. It is a post-freeze discovery review layer.

## Frozen universe before broad RT results

- 999 candidates from 402 genomes
- 24,020 linked predicted ORFs
- 17,655 unique amino-acid sequences
- 35,150 high-complexity repeat-group descriptions whose full span is outside predicted CDS
- 11,788 unique repeat cores

Architecture strata are fixed before the broader RT scan:

- N1: exact high-complexity noncoding repeat group overlaps the mathematical candidate
- N2: approximate high-complexity noncoding repeat group overlaps the mathematical candidate
- N3: auxiliary RT-like profile context
- N4: remaining repeat-context candidates

Every one of the 999 candidates is analyzed; strata control review order only.

## Broad RT branch

Primary screen uses the published **myRT RVT-All.hmm** class-specific RT profile database with explicit `hmmscan -E 1e-4 --domE 1e-4`. myRT reports 45 HMM models covering 41 RT classes. Primary hits are then searched against the myRT `cdd-pfamA.hmm` collection. A fixed text rule marks RT-related refinement domains; raw domain tables are retained.

A separate `YxDD` motif-only branch is descriptive and cannot establish RT identity.

## Registered-reference proximity

Three previously registered ART RTs, three registered ART partner references, and three registered array cores are searched only to flag obvious proximity to these references. This is **not** a complete comparison to Anthropic's 95-member ART catalogue and cannot establish novelty by itself.

Heuristic proximity flags are intentionally conservative and descriptive:

- RT: >=60% identity and >=70% coverage on both sequences
- partner: >=40% identity and >=70% coverage on both sequences
- array core: >=80% identity and >=80% coverage on both sequences

The final `T0_HIGH_NOVELTY_REVIEW_NOT_CONFIRMED` tier means only: broad RT support + noncoding array architecture + no close component among these three registered reference sets. It is **not** `new ART confirmed` and is not yet wet-lab-ready.

## Run

Recommended isolated environment:

```bash
conda create -n art_novelty -y --override-channels --strict-channel-priority \
  -c conda-forge -c bioconda \
  python=3.11 pandas biopython=1.86 hmmer=3.4 blast=2.16.0
conda activate art_novelty
```

Then:

```bash
./run_selftest.sh
./scripts/prepare_references.sh
./run_stage05d1.sh 8
cat outputs/GATE_STAGE05D1.json
cat outputs/SUMMARY_STAGE05D1.json
./collect_results.sh
```

No UG27 files are read or written.
