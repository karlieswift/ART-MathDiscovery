# ART MathDiscovery Stage05D4-R5 — Boltz-2 full-3D P1

This package freezes the first full-coordinate structure-computation tier after Stage05D4-R4.

## Scientific boundary
P1 is **not** a biological candidate filter. It is a compute-priority tier inside the R4 unresolved set. All 244 R4 AF/Swiss-Prot 3Di no-hit proteins, all 350 R3 Tier-A proteins, and all 602 DNA-first N1/N2 candidates remain scientifically retained.

P1 is prespecified before Boltz results: R4 no-hit; 220–500 aa; minimum linked ORF-to-mathematical-candidate distance <=250 bp; linked to >=2 DNA candidate loci. Expected P1=43. P2 uses the same rules but singleton candidate linkage; expected P2=40 and is not run in this package.

Primary prediction uses Boltz-2 with remote MSA server, fixed seed 20261007, default 3 recycling / 200 sampling steps, 1 diffusion sample. There is **no automatic single-sequence fallback**. If MSA networking fails, that is a technical failure; an `msa: empty` branch would be a separately labeled sensitivity, never silently mixed into primary results.

After Boltz-2, full-coordinate Foldseek searches predicted structures against the already frozen Stage05D4-R4 AlphaFold/Swiss-Prot structure DB. Neither Boltz confidence nor a Foldseek hit assigns RT/ART identity. No-hit is not negative.

## Run
```bash
conda activate art_structure
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
R3="$P/ART_MATH_DISCOVERY_STAGE05D4R3_COMPACT_STRUCTURE_READY_v1_0_1_20261007"
R4="$P/ART_MATH_DISCOVERY_STAGE05D4R4_PROSTT5_FOLDSEEK_REMOTE_STRUCTURE_v1_0_0_20261007"
./run_selftest.sh "$R3" "$R4"
./freeze_stage05d4r5.sh "$R3" "$R4"
./run_smoke_one.sh
./run_stage05d4r5_p1.sh 8

# Foldseek is in art_novelty, so either install foldseek in art_structure or run this step after `conda activate art_novelty`.
conda activate art_novelty
./run_full3d_foldseek_and_summarize.sh "$R4" 32
./collect_stage05d4r5.sh
```
