#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then echo "usage: $0 E3E_RESULTS_TAR PROJECT_ROOT" >&2; exit 2; fi
E3E=$1; P=$2
python scripts/stage05e3f.py --e3e "$E3E" --project-root "$P" --outdir outputs --workdir work
