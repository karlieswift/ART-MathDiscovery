# Stage04B-R1D2: RT / ORF / exact-core repeat evidence

Version: v1.0.0, 2026-10-05. This is a new post-freeze biological-evidence
stage, not an accession re-download, mathematical rescore, or replacement
for the audited R1D1 official proteins. All 383 candidates are retained.

## What is already inside this package

* The exact frozen 383-candidate registry and all original 25 fields.
* All 174 source-verified genome FASTAs, derived from the frozen Datasets cache.
* The exact 9,576-record official neighborhood protein FASTA and registry.
* A complete 20,804-node official feature ledger for CDS overlap checks.
* A completed **official-CDS-only** repeat/context preview for all 383
  neighborhoods, under the frozen method policy, with independent validation.
* Executable code, tests, locks and collection script.

No upstream NCBI download or project-root lookup is required. Original cache
archives, prior Stage04A outputs and UG27 runtime are not modified or accessed.

**The HMMER and Prodigal candidate scans have NOT been run by the package
builder. The full RT/ORF evidence gate is produced only by the server run.**
The three Pfam reference HMMs are fetched and byte-frozen before that run;
they are not fabricated or pre-included as unverified models.

## Install into an isolated environment

Do not install into or change `bioh200`, `drhu` or `blast_tools`.

```bash
conda create -n art_rt_evidence -y \
  --override-channels --strict-channel-priority \
  -c conda-forge -c bioconda \
  python=3.11 hmmer=3.4 prodigal=2.6.3

conda activate art_rt_evidence
python3 --version
hmmsearch -h | head -4
prodigal -v
```

The Python implementation uses the standard library only. It does not
require Biopython, NumPy, pandas, pyhmmer or pyrodigal at runtime.

## Run

Put the tar.gz and its sha256 file in your ART_MathDiscovery_v1 root.

```bash
cd /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1

sha256sum -c \
  ART_MATH_DISCOVERY_STAGE04B_R1D2_RT_ORF_ARRAY_EVIDENCE_v1_0_0_20261005.tar.gz.sha256 &&

tar -xzf \
  ART_MATH_DISCOVERY_STAGE04B_R1D2_RT_ORF_ARRAY_EVIDENCE_v1_0_0_20261005.tar.gz &&

cd ART_MATH_DISCOVERY_STAGE04B_R1D2_RT_ORF_ARRAY_EVIDENCE_v1_0_0_20261005 &&

./run_selftest.sh &&
./prepare_references.sh &&
./run_stage04b_r1d2.sh 4
```

Use a single run directory and a single process. `flock` prevents duplicate
writers. Four CPUs are a conservative operational choice while UG27 runs.
No claim about wall-clock time or maximum RAM is made.

The reference step downloads only PF00078, PF07727 and PF13456 through the
InterPro/Pfam API. It validates model identity, amino-acid alphabet,
accession.version, gathering thresholds and completeness, then freezes
model SHA256, source and a combined panel. A failed fetch stops before any
candidate score is computed. It never silently swaps models or relaxes cutoffs.

For an already available, authoritative local Pfam-A model library, the
same fixed three-model panel can instead be extracted before reference freeze:

```bash
./prepare_references.sh --pfam-hmm /absolute/path/Pfam-A.hmm.gz
```

Do not change reference libraries after `refs/REFERENCE_LOCK.json` exists.
The local filename in that optional command must be replaced with a real
existing path. Downloading an entire Pfam library is not the default workflow.

## Stages actually executed by the run command

1. Check all input hashes and 383/174/9576/9330 cardinalities.
2. Check HMMER 3.4 / Prodigal 2.6.3 and freeze tools, inputs, code, policy
   and reference versions in `outputs/RUN_LOCK.json`.
3. Exercise each real HMM using its consensus as a technical smoke test.
   This is NOT a biological sensitivity or specificity validation.
4. Predict coding ORFs on every one of the 174 complete genomes with
   Prodigal meta mode, and verify predicted CDS nucleotides against the
   exact genome coordinates and predicted coding direction.
5. Keep the genome-only rescue branch and the 173-genome supplemental ORF
   branch explicitly separate from official protein evidence.
6. Select predicted ORFs by overlap with the already fixed candidate
   neighborhoods. Freeze one combined exact-AA-deduplicated target universe;
   retain every official/predicted record and every candidate link.
7. Run the fixed Pfam panel using both GA and the declared sensitivity tier.
8. Compute exact-core direct-repeat descriptors for all 383 neighborhoods,
   independently of whether any RT profile hit occurred. Check the entire
   repeat span and individual copies against official and predicted CDS masks.
9. Integrate results back into all 383 candidates without changing any
   original field or dropping any candidate. Report RT-repeat distances
   and neighboring ORFs, without claiming functional association or novelty.

See `docs/METHODS_AND_BOUNDARIES.md` for the exact scientific scope.

## Outputs

```text
outputs/GATE_STAGE04B_R1D2.json
outputs/SUMMARY_STAGE04B_R1D2.json
outputs/ALL_383_RT_ORF_REPEAT_EVIDENCE.tsv
outputs/RT_CORE_GA_CANDIDATE_VIEW.tsv
outputs/ALL_PROTEIN_RT_ASSESSMENTS.tsv
outputs/RT_ALL_REPORTED_DOMAINS.tsv
outputs/KNOWN_RT_TEXT_CASE_SCAN_STATUS.tsv
outputs/PRODIGAL_ALL174_GENE_REGISTRY.tsv
outputs/PRODIGAL_NEIGHBORHOOD_REGISTRY.tsv
outputs/EXACT_CORE_REPEAT_GROUPS.tsv
outputs/RECURRENCE_CODING_CONTEXT.tsv
outputs/RT_REPEAT_PAIR_CONTEXT.tsv
outputs/RT_NEIGHBOR_ORF_INVENTORY.tsv
outputs/TARGET_UNIVERSE_LOCK.json
outputs/hmmer/
```

`RT_CORE_GA_CANDIDATE_VIEW.tsv` is a convenience view, not a survival cutoff.
The authoritative candidate denominator remains 383 in `ALL_383_...tsv`.
The output fields distinguish official RT hits, genome-only ORF rescue
hits and supplemental ORF hits. No domain hit is called an active RT or ART.

A successful full computation writes:

```text
PASS_ALL383_RT_ORF_EXACT_REPEAT_EVIDENCE_COMPUTED_WITH_SCOPE_LIMITS
```

An input/tool/parse/execution problem writes a HOLD gate and exits nonzero.
A HOLD, missing annotation, missing profile hit or missing exact repeat is
never converted into a negative biological call.

## Resume and collect

Verified Prodigal accession outputs and completed HMMER scans are reused
only with matching input/code/tool/reference locks and output hashes.
Rerun `./run_stage04b_r1d2.sh 4` in the SAME directory after interruption.
Do not edit its inputs, policy, references or scripts between runs.

After completion, whether PASS or HOLD:

```bash
cat outputs/GATE_STAGE04B_R1D2.json
./collect_results.sh
```

Send the timestamped `*_RESULTS_*.tar.gz` and `.sha256`. The collector includes
raw ORF predictions and HMMER outputs but does not duplicate the bundled
input genomes, upstream 174-cache or any UG27 arrays.

## Already computed preview is not the final RT/ORF result

`evidence/official_only_repeat_preview/` contains the actual local repeat
preview, 383 candidate-level context rows and 851 neighborhood-level repeat
records. There are 648 unique genomic repeat-group descriptions after exact
coordinate/core deduplication; neighborhoods overlap. No RT/ORF result is
inferred from these counts. The reference and server-run output directories
remain unpopulated in the delivery archive.
