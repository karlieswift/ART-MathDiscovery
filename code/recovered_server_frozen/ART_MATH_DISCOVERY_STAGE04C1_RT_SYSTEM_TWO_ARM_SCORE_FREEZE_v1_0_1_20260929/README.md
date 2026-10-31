# Stage04C1 RT-system two-arm score freeze v1.0.1

This supersedes v1.0.0 before scientific scoring. v1.0.0 hard-coded the obsolete Stage04C0-R1 manifest/gate filenames even though R2/R3 sequence-source corrections were required.

v1.0.1 requires a PASS Stage04C0-R3 sequence freeze and scans only:
- SEN2_14028S
- EC78_ECO7
- DRT2_KPN

Scoring remains the exact frozen two-arm logic: ACNV family + generic recurrence, q<=0.005, strict-overlap union, no Top-K. No RT/ncRNA/partner coordinates are read before score freeze.
