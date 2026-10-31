# Stage05D4 — N1/N2 divergent-RT remote-homology preparation

Purpose: freeze and prepare the Stage05D1 N1/N2 candidate-linked protein universe for divergent RT follow-up without changing any DNA score/rank/candidate.

Primary scientific contract:
- Source candidate universe is Stage05D1 only.
- N1 + N2 candidate count must equal 602 (262 + 340).
- DNA scores/ranks are never recomputed.
- Candidate filtering is forbidden; protein-level labels are post-discovery interpretation only.
- Formal broad-RVT positives are marked, not used to delete candidate loci.
- Exact-AA redundancy is collapsed only for compute efficiency; fanout is preserved.
- No claim of RT/ART identity is made in this stage.

Run:

```bash
conda activate art_novelty
D1=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05D1_NOVELTY_FIRST_MINING_v1_0_0_20261006
./run_selftest.sh "$D1"
./run_stage05d4.sh "$D1"
./collect_stage05d4.sh
```
