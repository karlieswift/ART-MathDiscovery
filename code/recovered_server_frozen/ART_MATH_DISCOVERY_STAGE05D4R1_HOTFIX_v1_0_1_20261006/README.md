# Stage05D4-R1 v1.0.1 hotfix

Pure implementation hotfix. MMseqs FASTA/cluster identifiers are `AA_<sha256>` while the Stage05D4 protein-universe table stores bare `aa_sha256`. v1.0.0 joined the two namespaces directly and raised `KeyError` after both NR90 and FAM50 clustering had completed.

v1.0.1 strips only the `AA_` prefix at the table join boundary. It does not rerun or alter MMseqs clustering, thresholds, the 602-candidate universe, DNA scores/ranks, or any scientific selection rule.

Run from the existing v1.0.0 Stage05D4-R1 directory after extracting this hotfix elsewhere:

`/path/to/hotfix/apply_hotfix_and_finish.sh "$D4" "$PWD"`
