#!/usr/bin/env python3
from pathlib import Path
from common import assert_manifest
root=Path(__file__).resolve().parents[1]
p=root/'PACKAGE_MANIFEST.tsv'
if not p.is_file():raise SystemExit('Package manifest missing')
assert_manifest(root,p)
print('[PASS] delivery package file hashes')
