#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 2 ]; then echo "usage: $0 F1_DIR F2B_DIR" >&2; exit 2; fi
mkdir -p outputs inputs/boltz2 work logs
python scripts/prepare_stage05f2c.py "$1" "$2" | tee logs/prepare.log
