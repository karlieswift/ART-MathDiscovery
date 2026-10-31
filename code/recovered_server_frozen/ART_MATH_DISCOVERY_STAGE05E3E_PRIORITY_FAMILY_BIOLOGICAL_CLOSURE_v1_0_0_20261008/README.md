# ART Mathematical Genome Discovery — Stage05E3E v1.0.0

Priority-family biological closure for ARCH250_0010 / ARCH250_0013 / ARCH250_0033.

The stage uses the already-frozen Stage05E3D v1.0.5 sequences and the already-existing Stage05C1 full 919-genome predicted proteome. It runs **phmmer only for local homology expansion**, then evaluates exact frozen repeat-core colocalization and topology around those homologs. If the existing `RT_LOCI_ALL919.tsv` can be resolved, it also adds RT proximity descriptively. It never changes E1/E2/E3 membership or DNA scores/ranks.

Recommended environment: `art_rt_evidence`, because the project has already used HMMER 3.4 there.

Run:

```bash
conda activate art_rt_evidence
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"

sha256sum -c ART_MATH_DISCOVERY_STAGE05E3E_PRIORITY_FAMILY_BIOLOGICAL_CLOSURE_v1_0_0_20261008.tar.gz.sha256

tar -xzf ART_MATH_DISCOVERY_STAGE05E3E_PRIORITY_FAMILY_BIOLOGICAL_CLOSURE_v1_0_0_20261008.tar.gz

E3E="$P/ART_MATH_DISCOVERY_STAGE05E3E_PRIORITY_FAMILY_BIOLOGICAL_CLOSURE_v1_0_0_20261008"
E3D="$P/ART_STAGE05E3D_PREDICTED_WITNESS_HASH_RESCUE_RESULTS_20261008_131220.tar.gz"

cd "$E3E"
./tests/run_selftest.sh
./run_stage05e3e.sh "$E3D" "$P" 16

cat outputs/STAGE05E3E_GATE.json
column -t -s $'\t' outputs/STAGE05E3E_FAMILY_CLOSURE_SUMMARY.tsv
column -t -s $'\t' outputs/STAGE05E3E_EXACT_REPEAT_COLOCALIZATION_SUMMARY.tsv
column -t -s $'\t' outputs/STAGE05E3E_RT_CONTEXT_SUMMARY.tsv
cat outputs/STAGE05E3E_MANUSCRIPT_NOTES.md
./collect_stage05e3e.sh
```

Return the result tar.gz and `.sha256`.

Important: phmmer hits outside the frozen E2A/E3D witness set are **local homologs**, not newly admitted family members. Exact 12–20-nt repeat-core matches near a homolog are architecture-support observations, not proof of function or global uniqueness.
