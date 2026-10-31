# ART MathDiscovery Stage05D4-R4

Fast orthogonal remote-structure triage of the frozen 350 Stage05D4-R3 Tier-A proteins using Foldseek + ProstT5 against AlphaFold/Swiss-Prot.

This stage does **not** recompute DNA scores/ranks, does not filter the 602 N1/N2 candidate universe, and does not treat either a structural hit or a no-hit as an ART identity call.

## Run

1. Install Foldseek if absent (Bioconda is sufficient for CPU ProstT5 mode):
   `conda install -n art_novelty -c conda-forge -c bioconda foldseek -y`
2. `./run_selftest.sh "$R3"`
3. `./prepare_foldseek_refs.sh 32`
4. `./run_stage05d4r4.sh "$R3" 32 0`
5. `./collect_stage05d4r4.sh`

Use GPU=1 only if the installed Foldseek binary has CUDA/ProstT5 GPU support. CPU is the safe default.
