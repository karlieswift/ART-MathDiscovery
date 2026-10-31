# Stage05B: full 919-new-genome frozen two-arm discovery scan

This package advances **new ART-like system discovery for later wet-lab follow-up**.
It scores DNA; it does not claim an RT, an RNA array, a partner function, or a novel
ART identity before those layers have been examined.

## Frozen source cohort

- Metadata selected: 980 (250 / 250 / 480 across three length bins).
- Sequence-eligible, non-exact-duplicate genomes: **919** (237 / 233 / 449).
- 171,397,507 bp; 21,062,200 oriented windows.
- Exclusion ledger: 46 non-ACGT; 13 exact/RC matches to prior240; 2 exact/RC matches
  within the new selection. No replacement after these exclusions.
- Original Stage05A manifest SHA256:
  `5eec6d7ee00584c6266b8a953b6112540148b949341e400099b285ce30d235f0`.

The returned Stage05A package contains metadata, manifests and ledgers, not the
919 raw FASTAs. Therefore this package **rechecks all 919 raw source files on your
server before starting scoring**. No redownload, base repair, sequence trimming,
or replacement accession is permitted.

## Engine choice and reproducibility

The scientific engine is the already-executed **Stage04C1 clean two-arm v1.0.2**
from 2026-09-29. Its entrypoint, common functions, ACNV extractor and policy file
are copied byte-for-byte. This is not a new detector trained or tuned on 919
sequences. ACNV and generic recurrence are computed from DNA; recurrence does
not consume ACNV vectors. All six window lengths, quarter-window strides,
strands, family calibration, q=0.005, strict-overlap union and priority rules
remain those of that engine.

This package does **not** invoke the running UG27 job, does not depend on its
final gate, and does not read UG27 locus truth. It does not change the prior240
boards or the prior383 registry. The pilot's historical region IDs/ranks are not
recomputed or silently redefined as the new clean-two-arm board.

Internal filenames starting `C1_` and the original `control_id` field are retained
for exact source-lineage preservation. The current scientific experiment is
Stage05B new-source discovery; those names do not make the 919 genomes controls.

## Run (Linux, Python >=3.9, NumPy and pandas)

Use the existing `bioh200` environment on as-4125. No BLAST, HMMER, Prodigal,
Datasets CLI, GPU, network call, or new dependency installation is needed.

```bash
conda activate bioh200
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"
sha256sum -c ART_MATH_DISCOVERY_STAGE05B_NEW919_FROZEN_TWO_ARM_SCAN_v1_0_0_20261005.tar.gz.sha256
tar -xzf ART_MATH_DISCOVERY_STAGE05B_NEW919_FROZEN_TWO_ARM_SCAN_v1_0_0_20261005.tar.gz
cd ART_MATH_DISCOVERY_STAGE05B_NEW919_FROZEN_TWO_ARM_SCAN_v1_0_0_20261005
./run_selftest.sh
./run_stage05b.sh "$P/ART_MATH_DISCOVERY_STAGE05A_NEW_SOURCE_DISCOVERY_FREEZE_v1_0_0_20261005" 2 8
```

The one command runs: (1) local source-byte preflight; (2) two historical small
whole-genome replays and exact checks; (3) the entire frozen 919-genome queue;
(4) all-genome and all-region final score freeze. Historical failures stop the
new scan automatically. There is no required additional chat round to authorize
starting the entire cohort after successful preflight/validation.

For a plan-only check: `./run_preflight.sh STAGE05A_ROOT 2 8`.
For a historical-only check: `./run_historical_validation.sh 8`.
For a read-only progress view: `./status.sh`.

## Resource and restart scope

The 761.44 GB number printed by Stage05A was the **cumulative** legacy feature
array footprint across all genomes, not peak memory/disk. The already-frozen
clean-two-arm source uses one 5,468-feature float32 ACNV array per genome.
Largest such array: **1,483,927,712 bytes**; two concurrent genomes are bounded
by **2,967,855,424 bytes of these arrays**. Total generated ACNV array bytes over
the full queue are 460,672,438,400; arrays are removed after each successful
genome. Score-table files and robust-score temporaries require additional disk
and RAM. These array figures are not process-RAM upper bounds.

Default: 2 simultaneous genomes x 8 feature worker processes. BLAS inner
threading is capped at 1. The preflight checks CPU affinity, a conservative
8 GiB available RAM per genome worker plus 4 GiB headroom, and 20 GiB free disk.
It cannot reserve CPU for unrelated jobs; keep UG27's existing resources intact.
Workers/threads are operational settings and may be reduced on restart without
changing any genome, window, score definition or queue membership.

Genomes are ordered by length descending, then accession, only for scheduling.
This is not a biological rank, a large-genome-only subset, or adaptive selection.
All 919 must complete before the cohort gate passes.

**Restart:** rerun the same command. Completed genome boards are verified and
reused; an interrupted genome's feature extraction starts again from its input.
This release does not claim within-genome stratum checkpoints. A file lock
prevents another writer or collector in the same package directory. Different
Python/NumPy/pandas runtimes cannot silently be mixed into the same run.

## Historical validation

The two included references are `NC_002666.1` and `NC_028671.2`. Their raw DNA was
reconstructed from frozen cache strings using the historical header/line wrap;
the entire FASTA file SHA256 matches the original Stage04A manifest.

The byte-frozen clean engine rescans these whole genomes. Every shared window
field (member scores, strand ties, family q, recurrence values) must exactly match
the original Stage04A score text. Clean-two-arm proposals and union boards must
exactly match those deterministically derived from that archived score table
using the frozen Stage04C1 merge. No numerical tolerance is used. Old pilot
region ranks are **not** the expected clean-two-arm rank reference.

## Output and handoff

- `outputs/PREFLIGHT_GATE.json`: all919 local FASTA-byte audit.
- `outputs/HISTORICAL_EQUIVALENCE_GATE.json`: exact historical checks.
- `outputs/STAGE05B_FINAL_GATE.json`: only PASS after all919 boards exist/verify.
- `outputs/STAGE05B_ALL_TWO_ARM_UNION_REGIONS.tsv`: all OR-supported regions.
- `outputs/STAGE05B_DUAL_ARM_CANDIDATES.tsv`: two-arm support view; no Top-K cut.
- `outputs/STAGE05B_PER_GENOME_BURDEN.tsv`: all919, including zero-candidate rows.
- `outputs/SCORE_FILE_MANIFEST.tsv`: full score-file commitments.
- `outputs/genomes/<accession>/`: original per-genome boards and completion lock.

After a completed or stopped run, `./collect_results.sh` creates a result TAR and
SHA256. By default it includes region/proposal tables, per-genome gates,
aggregate outputs, logs, code and provenance, but **omits raw genomes, numeric
matrices, and bulky window-score tables**. Omitted window files remain local;
their hashes are listed. `./collect_results.sh --with-window-scores` explicitly
includes the window-score bulk when required.

Archive and hash the completed score board before biological annotation.
Next: RT/ORF evidence, imperfect/structural array evidence calibrated separately,
CDS-versus-intergenic context, homologous-system dereplication, known-system
comparison, and candidate evidence cards for experimental prioritization.
No fixed count of wet-lab candidates is promised.
