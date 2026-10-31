# ART Mathematical Genome Discovery - Stage05E3B

Family-level lineage dereplication + official-context synthesis for all 12 frozen E2A exact-primary joint-signature families.

Primary input: `ART_STAGE05E3A_E2A_OFFICIAL_CONTEXT_LINEAGE_RESULTS_R1_20261007_172658.tar.gz` (SHA256 `31c70ee1c7769baec58da11dfcc40fc3ec6ab7c17d6608391f737bc404041bca`).

Run:

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P/ART_MATH_DISCOVERY_STAGE05E3B_FAMILY_DEREPLICATION_v1_0_0_20261007"
./tests/run_selftest.sh
./run_stage05e3b.sh "$P/ART_STAGE05E3A_E2A_OFFICIAL_CONTEXT_LINEAGE_RESULTS_R1_20261007_172658.tar.gz" "$P"
cat outputs/STAGE05E3B_SUMMARY.json
column -t -s $'\t' outputs/STAGE05E3B_REVIEW_QUEUE.tsv | head -20
./collect_stage05e3b.sh
```

Expected frozen regression:
- 12 families, 119 candidates, 61 genomes, 80 exact-primary witness candidates.
- 64 official-supported witness proteins; 63/64 hypothetical/uncharacterized; median 64 aa; 52/64 <=100 aa.
- lineage breadth: 4 broad-multi-lineage, 4 multi-lineage, 4 narrow-replicated.
- review strata: E3B1=1, E3B2=2, E3B3=5, E3B4=2, E3B5=2.
- strict RT keyword families=0; any known-system-context family=1 (distant Cas4/CRISPR context in ARCH250_0036).

Gate: `PASS_E2A_FAMILY_LINEAGE_DEREPLICATION_AND_CONTEXT_SYNTHESIS_FROZEN`.
