# Stage04C0-R3 sequence-length correction

This hotfix supersedes only the Stage04C0-R2 sequence-length expectation for DRT2_KPN.

The R2 source mapping `NC_012731.1 -> AP006725.1` is correct. The failed R2 expected length `5,472,672 bp` was the sum of the NTUH-K2044 chromosome (`5,248,520 bp`) and pK2044 plasmid (`224,152 bp`). The frozen DRT2 control record is the chromosome `NC_012731.1/AP006725.1`, so the correct expected sequence length is `5,248,520 bp`.

R3 performs **no network download** and **no ACNV/recurrence scoring**. It copies/reflinks the already downloaded R2 records into its own directory, revalidates all four against exact chromosome lengths and writes a new frozen manifest.
