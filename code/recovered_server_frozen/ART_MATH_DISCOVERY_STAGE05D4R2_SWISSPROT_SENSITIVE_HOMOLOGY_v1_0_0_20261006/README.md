# Stage05D4-R2: Swiss-Prot sensitive homology triage

Input: frozen Stage05D4-R1 NR90 representatives (9,245 queries) derived from 602 N1/N2 DNA-first candidates.

Run:
```bash
./run_selftest.sh
./run_stage05d4r2.sh /path/to/ART_MATH_DISCOVERY_STAGE05D4R1_SEQUENCE_FAMILY_COMPRESSION_v1_0_0_20261006 32
```
The runner downloads the current UniProtKB/Swiss-Prot FASTA once, freezes SHA256, performs sensitive MMseqs2 search (`-s 7.5`, E<=1e-3, up to 25 hits/query), and writes query/candidate-level triage outputs. No-hit is not negative. No DNA score/rank/candidate is changed.
