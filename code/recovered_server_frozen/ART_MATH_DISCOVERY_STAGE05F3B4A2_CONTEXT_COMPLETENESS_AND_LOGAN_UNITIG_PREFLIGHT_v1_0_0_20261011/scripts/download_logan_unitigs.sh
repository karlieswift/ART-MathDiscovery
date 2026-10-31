#!/usr/bin/env bash
set -euo pipefail
ACC_FILE=$1
DEST=$2
WORKERS=${3:-2}
mkdir -p "$DEST"
command -v wget >/dev/null 2>&1 || { echo '[FATAL] wget missing' >&2; exit 2; }
export DEST
cat "$ACC_FILE" | awk 'NF' | xargs -P "$WORKERS" -I{} bash -c '
  a="$1"; f="$DEST/$a.unitigs.fa.zst"; u="https://s3.amazonaws.com/logan-pub/u/$a/$a.unitigs.fa.zst";
  echo "[DOWNLOAD] $a";
  wget -c "$u" -O "$f";
  test -s "$f";
  echo "[PASS] $a $(stat -c %s "$f") bytes";
' _ {}
