# ART MathDiscovery Stage05E3H v1.0.0

Accelerated global-novelty + short-protein functional-closure package for ARCH250_0010, ARCH250_0013 and ARCH250_0033.

Primary run uses only frozen E3D/E3F/E3G result archives and existing project artifacts. It does not need GPU.

Optional external resources can be supplied through environment variables:
- `E3H_PFAM_HMM=/path/to/Pfam-A.hmm` (pressed database preferred)
- `E3H_PHROG_FASTA=/path/to/phrog_proteins.faa`
- `E3H_VOG_FASTA=/path/to/vog_proteins.faa`
- `E3H_REMOTE_BLAST=1` to attempt NCBI remote BLASTP for the three representatives (off by default; failure is nonfatal)

The stage always emits a precise resource/status ledger so missing external databases remain explicit rather than being interpreted as negative evidence.
