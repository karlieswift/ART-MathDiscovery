#!/usr/bin/env bash
set -euo pipefail
J2_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ $# -lt 1 ]]; then
  printf '%s\n' 'Usage: bash run_j2.sh {preflight|admit|freeze|score|evaluate|historical|test|fetch-rvt} [arguments]'
  exit 2
fi
J2_COMMAND="$1"
shift
case "$J2_COMMAND" in
  historical) exec python "$J2_HOME/scripts/historical_smoke.py" "$@" ;;
  test) exec python "$J2_HOME/scripts/test_guards.py" "$@" ;;
  fetch-rvt) exec python "$J2_HOME/scripts/fetch_frozen_rvt.py" "$@" ;;
  preflight|admit|freeze|score|evaluate) exec python "$J2_HOME/scripts/j2_pipeline.py" "$J2_COMMAND" "$@" ;;
  *) printf '%s\n' "Unknown command: $J2_COMMAND" >&2; exit 2 ;;
esac
