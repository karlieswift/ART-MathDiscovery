# Stage04B9 — strict mirror + ARCH5 structure closure

Two independent branches are packaged together but must be reported separately.

## A. Strict mirror

54G raw/frozen delta -> 88,487 formal -> 85,503 exact-AA -> native Stage03 26,396 operational families -> exact historical F3B3A 4x8 rule -> 32 -> same Boltz2 -> same clean-RT Foldseek -> corrected set-level classification.

No Stage04A combined OLD/DELTA family filter is used. The package will HOLD rather than silently substitute the later B5 reclustering if the native Stage03 `F3V2_FAMILY_cluster.tsv` is unavailable.

## B. ARCH5

The five Stage04B8R1 unitig-ABSENT candidates proven complete stop-to-stop on their original source contigs are structurally tested under the same frozen protocol. They remain a post-B8 technical-completeness extension and are not retroactively added to the original B7 primary 12.

## Run

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P/ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B9_STRICT_MIRROR_AND_ARCH5_STRUCTURE_CLOSURE_v1_0_0_20261021"
./run_stage04b9_prepare.sh "$P"
```

Then run the short ARCH5 branch first:

```bash
conda activate art_structure
./run_stage04b9_boltz_arch5.sh "$P" 0,1,2,3
```

Then the strict mirror (existing exact-AA B3 structures are reused only when unambiguous; only missing structures are predicted):

```bash
./run_stage04b9_boltz_strict_mirror.sh "$P" 0,1,2,3
```

Both runners are resume-safe: existing CIFs are skipped. If a GPU is occupied and a query OOMs, rerun with a freer GPU list; do not change scientific parameters.

Finally:

```bash
conda activate art_novelty
./run_stage04b9_foldseek_finalize.sh "$P" 32
```

Key outputs:
- `outputs/STAGE04B9A_PREPARE_GATE.json`
- `outputs/STRICT_MIRROR_SHORTLIST_32.tsv`
- `outputs/STRICT_MIRROR_VS_DELTAONLY32_EXACT_AA_OVERLAP.tsv`
- `outputs/STRICT_MIRROR_STRUCTURE_GATE.json`
- `outputs/ARCH5_STRUCTURE_GATE.json`
- `outputs/ARCHITECTURE_FIRST_DELTAONLY17_COMPLETE_STRUCTURAL_CLOSURE.tsv`
- `outputs/STAGE04B9_CLOSEOUT_GATE.json`
- `docs/SUPPLEMENT_ROUTE_TABLE_TEMPLATE.md`
