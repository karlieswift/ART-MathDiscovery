#!/usr/bin/env bash
set -euo pipefail
OLD=${1:?path to completed v1.0.0 package directory required}
HERE=$(cd "$(dirname "$0")" && pwd)
[ -d "$OLD/work/accession_results" ] || { echo "HOLD missing $OLD/work/accession_results"; exit 2; }
mkdir -p "$HERE/work"
# Copy the completed scientific scan products. These are small relative to raw 54G data.
cp -a "$OLD/work/." "$HERE/work/"
echo "PASS_STAGE04B5_V101_IMPORTED_V100_WORK"
find "$HERE/work/accession_results" -mindepth 1 -maxdepth 1 -type d | wc -l | awk '{print "imported_accession_result_dirs=" $1}'
