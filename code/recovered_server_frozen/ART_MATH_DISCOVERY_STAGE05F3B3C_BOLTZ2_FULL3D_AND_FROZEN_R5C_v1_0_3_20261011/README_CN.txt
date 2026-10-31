F3B3C v1.0.3: CRLF/Foldseek handoff implementation hotfix.

Root cause in v1.0.2:
Python csv.DictWriter used CRLF. Bash TSV parsing retained the trailing carriage return in cif_path, so Foldseek received a non-existent path ending in \r. set -e/pipefail then exited at the first Foldseek call; no .m8 result was produced.

This package:
- reuses v1.0.2 Boltz outputs/status; no Boltz recomputation;
- normalizes the Foldseek job handoff to LF;
- keeps frozen R5C DB and thresholds unchanged;
- records Foldseek technical failures explicitly so they can never be misclassified as NO_HIT.
