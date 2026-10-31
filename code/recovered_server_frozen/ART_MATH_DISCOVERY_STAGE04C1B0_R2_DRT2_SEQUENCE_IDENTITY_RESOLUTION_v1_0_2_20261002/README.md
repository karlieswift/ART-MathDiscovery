# Stage04C1B0-R2 DRT2 sequence-identity truth-resolution hotfix

This hotfix supersedes the failed Stage04C1B0-R1 v1.0.1 truth-resolution attempt only. It does not read, recompute, mutate, or rank any Stage04C1 score.

R1 successfully replaced the fragile Ec78 product-label heuristic with the frozen sequence-level rule (unique 311-aa CDS with exact YADD), but then failed at DRT2 because the AP006725.1 INSDC source record does not need to expose the modern RefSeq protein accession `WP_012737279.1` or a modern `reverse transcriptase` product string. The exact frozen DRT2 CDS itself was present at 598753..600030(-).

R2 therefore verifies DRT2 sequence-first: exact frozen CDS coordinates and strand, exactly 425 translated amino acids, exact catalytic `YCDD` at aa 267-270 (D269/D270), and the frozen AP006725.1 source-translation SHA256. It does not require source-record product text or a `WP_` protein_id. Ec78 and Sen2 rules are unchanged from R1.

The already-downloaded GenBank flatfiles from the failed v1.0.0 directory are reused; no network fetch and no score access occur.
