# Stage05E3G Priority-family orthogonal evidence recovery

Input: Stage05E3F result tar + project root.

Run:

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
E3G="$P/ART_MATH_DISCOVERY_STAGE05E3G_PRIORITY_FAMILY_ORTHOGONAL_EVIDENCE_RECOVERY_v1_0_0_20261008"
E3F="$P/ART_STAGE05E3F_PRIORITY_ARCHITECTURE_EXPANSION_RT_CONTEXT_RESULTS_20261008_144517.tar.gz"
cd "$E3G"
./tests/run_selftest.sh
./run_stage05e3g.sh "$E3F" "$P"
cat outputs/STAGE05E3G_GATE.json
cat outputs/STAGE05E3G_FAMILY_REPRESENTATIVES.tsv
cat outputs/STAGE05E3G_STAGE_EVIDENCE_SUMMARY.tsv
cat outputs/STAGE05E3G_SHORT_PROTEIN_EXTENSION_PLAN.tsv
cat outputs/STAGE05E3G_MANUSCRIPT_NOTES.md
./collect_stage05e3g.sh
```

The stage is deliberately read-only with respect to earlier evidence. It does not re-run HMMER, MMseqs, Foldseek, Boltz, DNA scoring, or candidate ranking.
