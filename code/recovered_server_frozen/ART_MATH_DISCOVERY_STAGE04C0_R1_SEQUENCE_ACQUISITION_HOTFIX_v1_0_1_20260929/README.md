# Stage04C0-R1 sequence acquisition hotfix

Purpose: repair the Stage04C0 v1.0.0 source-acquisition bug in which NCBI E-utilities returned an HTML blocked-diagnostic page with HTTP 200 and the original validator misclassified it as sequence.

This hotfix does NOT mutate the original C0 directory. It retrieves INSDC nucleotide records from DDBJ getentry, validates FASTA/flatfile syntax and sequence length, freezes SHA256, and writes a C1-ready manifest.
