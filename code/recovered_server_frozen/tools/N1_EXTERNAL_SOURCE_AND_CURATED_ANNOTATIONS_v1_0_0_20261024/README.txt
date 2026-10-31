N1_EXTERNAL_SOURCE_AND_CURATED_ANNOTATIONS_v1_0_0_20261024

Purpose:
- public author/literature evidence screening, 12 source records
- five real exposed MarsHill nucleotide annotation verification
- prior exposure accession screening (923 historical accession versions)
- optional INPHARED2 metadata + GenBank RT-CDS candidate screening
- explicitly fail/hold if new independent native array+RT truth absent

Server:
  P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
  NAME=N1_EXTERNAL_SOURCE_AND_CURATED_ANNOTATIONS_v1_0_0_20261024
  mkdir -p "$P/tools"
  unzip -n "$P/$NAME.zip" -d "$P/tools"
  K="$P/tools/$NAME"
  python -m unittest discover -s "$K/tests" -v
  bash "$K/run_n1_external_screen.sh" "$P"

For independent external genome metadata:
  bash "$K/run_n1_external_screen.sh" "$P" --metadata /path/to/INPHARED2_selected.tsv
For author GenBank annotations:
  bash "$K/run_n1_external_screen.sh" "$P" --gbff-dir /path/to/new_author_gbff

To return already computed N3/N4/N5 results:
  bash "$K/pack_n345_replay_results.sh" "$P"

IMPORTANT: The external N1 primary benchmark remains HOLD. It requires
independent native genome array and RT truth, signed before looking at Arm M predictions.
