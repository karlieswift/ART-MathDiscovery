#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python -c 'import sys; sys.path.insert(0,"scripts"); from run_stage import verify_package; verify_package()'
python -m unittest discover -s tests -v
