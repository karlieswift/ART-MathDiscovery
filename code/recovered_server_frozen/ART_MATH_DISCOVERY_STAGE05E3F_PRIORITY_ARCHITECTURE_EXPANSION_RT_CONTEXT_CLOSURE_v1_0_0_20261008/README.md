# Stage05E3F — priority architecture expansion + RT-context closure

Purpose: read the already-PASS Stage05E3E result, separate protein-family expansion from exact repeat-linked architecture expansion for ARCH250_0010 / 0013 / 0033, and repair the previously unavailable Stage05C1 RT-context join without rescoring anything.

Key scientific expectations from the received E3E result, to be verified rather than assumed by E3F:
- ARCH250_0010: five exact E3D witnesses retain the 19-nt upstream repeat at 277 bp; one additional weak partial local homolog carries the exact repeat only at ~757 bp.
- ARCH250_0013: five non-witness local homologs add to the two E3D witnesses; all seven local homologs carry the exact 12-nt repeat within 500 bp on the upstream side.
- ARCH250_0033: a broad protein homolog family is present, but the exact 20-nt repeat is confined to the three E3D witness genomes in the E3E 2-kb search.

E3F treats these as descriptive within-project architecture patterns. It does not define global novelty or function.

Run:
```bash
conda activate art_rt_evidence
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P/ART_MATH_DISCOVERY_STAGE05E3F_PRIORITY_ARCHITECTURE_EXPANSION_RT_CONTEXT_CLOSURE_v1_0_0_20261008"
./tests/run_selftest.sh
./run_stage05e3f.sh \
  "$P/ART_STAGE05E3E_PRIORITY_FAMILY_BIOLOGICAL_CLOSURE_RESULTS_20261008_141923.tar.gz" \
  "$P"
cat outputs/STAGE05E3F_GATE.json
cat outputs/STAGE05E3F_FAMILY_EXPANSION_SUMMARY.tsv
cat outputs/STAGE05E3F_RT_CONTEXT_SUMMARY.tsv
cat outputs/STAGE05E3F_MANUSCRIPT_NOTES.md
./collect_stage05e3f.sh
```

If the E3E result tar is nested under the E3E package directory, basename-unique recovery is used. RT_LOCI_ALL919.tsv is first searched as a live project file; if absent, E3F searches the Stage05C1 result tar and extracts only that frozen table. The RT search is never rerun.
