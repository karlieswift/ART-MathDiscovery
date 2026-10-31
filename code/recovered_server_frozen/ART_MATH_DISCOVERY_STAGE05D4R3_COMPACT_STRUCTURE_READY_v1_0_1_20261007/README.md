# Stage05D4-R3 compact structure-ready freeze

This stage does **not** rescore DNA, rerank candidates, assign RT/ART identity, or remove any of the 602 frozen N1/N2 mathematical candidates.

It maps the frozen Stage05D4-R2 NR90 Swiss-Prot triage back through the Stage05D4-R1 clusters and Stage05D1 ORF/candidate links. It creates deterministic computational priority tiers for later structure/remote-homology work.

Primary structure tier (predeclared):
- Swiss-Prot label: `NO_SWISSPROT_HIT_NOT_NEGATIVE`
- representative protein length 180-700 aa
- minimum linked ORF-to-N1/N2-math-candidate distance <=1000 bp

Extended tier:
- same no-hit and length rule
- distance >1000 and <=2500 bp

No-hit is not interpreted as novelty or RT evidence. Tier A/B membership is only a compute-priority label.

Run:
```bash
./run_selftest.sh <D1> <D4> <R1> <R2>
./run_stage05d4r3.sh <D1> <D4> <R1> <R2>
./collect_stage05d4r3.sh
```

Implementation v1.0.1 normalizes the R2/architecture candidate-key namespace during table merge only; scientific tier rules are unchanged.
