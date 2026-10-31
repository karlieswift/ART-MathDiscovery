# Stage04B-R1A NCBI Datasets offline annotation cache

Purpose: retrieve exact versioned NCBI Datasets packages for the 174 sequences that carry the 383 frozen Stage04A-R1 two-arm candidate regions. This stage is retrieval/provenance only. It does not rescore, rerank, filter, or biologically interpret any candidate.

Canonical frozen registry SHA256:
`73007864fdef7dba662c9b90a02e2b915359c1fbd93a6061edeea35afde07349`

Expected accessions: 174.

Run on a host where `datasets` works (tested with NCBI Datasets CLI 18.38.0):

```bash
./run_selftest.sh
./run_stage04b_r1a_fetch.sh
./collect_stage04b_r1a.sh
```

If the gate PASSes, create the portable cache bundle:

```bash
./pack_cache.sh
```

Copy that cache bundle to the canonical analysis host; Stage04B-R1 biological annotation should consume this offline cache rather than contacting blocked EUtils/SViewer endpoints.
