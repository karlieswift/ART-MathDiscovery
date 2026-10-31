Stage05F3B3C v1.0.1 implementation-only hotfix
Date: 2026-10-11

Reason for v1.0.1:
- v1.0.0 stopped at preflight with: [FATAL] missing foldseek
- art_structure contains Boltz but does not contain Foldseek.
- Established project runtime uses Foldseek in art_novelty.
- No Boltz prediction, Foldseek result, structure tier, or scientific result was produced before the v1.0.0 failure.

v1.0.1 changes only runtime dispatch:
- Boltz-2 remains in art_structure.
- Foldseek is called through: conda run -n art_novelty foldseek
- Frozen 32-family handoff, Boltz settings, clean-RT DB, R5C thresholds and claim boundary are unchanged.

Run:
  conda activate art_structure
  export BOLTZ_CACHE=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache
  ./run_stage05f3b3c.sh 32 4

Do not truncate long proteins and do not substitute msa: empty on primary failures.
