# Stage05D4-R5 v1.0.7 downstream hotfix

This hotfix supersedes the failed v1.0.6 Foldseek invocation only.

Observed failure under Foldseek 10.941cd33:
`Format code tmscore does not exist.`

Cause: Foldseek supports `alntmscore` (and in current releases also query/target-normalized TM-score fields), not a generic `tmscore` output token.

v1.0.7 changes only the output-field spelling used by Foldseek:
`tmscore -> alntmscore`.

The existing Stage05D4-R5 summarizer still consumes column 11, so its legacy output columns `best_tmscore` and `best_rt_tmscore` specifically mean Foldseek `alntmscore` for v1.0.7.

The runner also rebuilds `work/p1_cif_flat_regular` from the 43 verified canonical regular Boltz CIFs before search, eliminating cwd/symlink dependencies.

No Boltz prediction is rerun. P1=43, target database, E-value 0.01, max-seqs 50, exact-tmscore=1, DNA scores/ranks, candidate universe and all discovery rules are unchanged.
