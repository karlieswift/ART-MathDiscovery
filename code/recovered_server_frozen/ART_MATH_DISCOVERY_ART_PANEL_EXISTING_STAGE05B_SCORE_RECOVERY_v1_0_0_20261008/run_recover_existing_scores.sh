#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo "usage: $0 <project_root> <ART_EXACT_SOURCE_REPLAY_INPUTS_*.tar.gz|dir> [Stage05B_results.tar.gz|dir]" >&2
  exit 2
fi
P="$1"; INP="$2"; S5="${3:-}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
rm -rf outputs work
mkdir -p outputs work logs
ARGS=(--project-root "$P" --exact-inputs "$INP" --outdir outputs --workdir work)
if [ -n "$S5" ]; then ARGS+=(--stage05b "$S5"); fi
python scripts/recover_stage05b_panel.py "${ARGS[@]}" 2>&1 | tee logs/recover.log
cat outputs/ART_PANEL_STAGEA1_GATE.json
