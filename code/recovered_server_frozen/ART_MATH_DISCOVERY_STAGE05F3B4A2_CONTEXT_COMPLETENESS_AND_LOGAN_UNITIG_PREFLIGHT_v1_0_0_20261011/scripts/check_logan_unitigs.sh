#!/usr/bin/env bash
set -euo pipefail
ACC_FILE=$1
OUT=$2
mkdir -p "$(dirname "$OUT")"
printf 'accession\turl\thttp_status\tcontent_length_bytes\tavailability\n' > "$OUT"
while IFS= read -r a; do
  [[ -n "$a" ]] || continue
  u="https://s3.amazonaws.com/logan-pub/u/$a/$a.unitigs.fa.zst"
  status=""; len=""
  if command -v curl >/dev/null 2>&1; then
    hdr=$(curl -L -sS -I --connect-timeout 15 --max-time 60 "$u" 2>/dev/null || true)
    status=$(printf '%s\n' "$hdr" | tr -d '\r' | awk '/^HTTP\//{s=$2} END{print s}')
    len=$(printf '%s\n' "$hdr" | tr -d '\r' | awk 'tolower($1)=="content-length:"{x=$2} END{print x}')
  elif command -v wget >/dev/null 2>&1; then
    hdr=$(wget --server-response --spider -T 60 "$u" 2>&1 || true)
    status=$(printf '%s\n' "$hdr" | awk '/HTTP\//{s=$2} END{print s}')
    len=$(printf '%s\n' "$hdr" | awk 'tolower($1)=="content-length:"{x=$2} END{print x}')
  else
    echo '[FATAL] need curl or wget' >&2; exit 2
  fi
  [[ -n "$status" ]] || status=NA
  [[ -n "$len" ]] || len=NA
  avail=NO
  [[ "$status" == "200" ]] && avail=YES
  printf '%s\t%s\t%s\t%s\t%s\n' "$a" "$u" "$status" "$len" "$avail" | tee -a "$OUT"
done < "$ACC_FILE"
