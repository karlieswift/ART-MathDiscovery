#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
python - "$HERE/scripts/graph_unique_path_recovery.py" <<'PY'
import importlib.util, sqlite3, sys
spec=importlib.util.spec_from_file_location('g',sys.argv[1]); g=importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
cx=sqlite3.connect(':memory:')
cx.executescript('CREATE TABLE nodes(id TEXT PRIMARY KEY, len INTEGER NOT NULL); CREATE TABLE edges(src_id TEXT, src_o TEXT, dst_id TEXT, dst_o TEXT);')
cx.executemany('INSERT INTO nodes VALUES (?,?)',[('u0',31),('u1',40),('u2',50),('b1',35),('b2',35),('p0',42)])
cx.executemany('INSERT INTO edges VALUES (?,?,?,?)',[('u0','+','u1','+'),('u1','+','u2','+'),('p0','+','u0','+')])
r=g.walk(cx,('u0','+'),'RIGHT',20)
assert r['status']=='FULL_UNIQUE_PATH' and r['recovered_bp']==30 and len(r['path'])==2, r
r=g.walk(cx,('u0','+'),'LEFT',10)
assert r['status']=='FULL_UNIQUE_PATH' and r['recovered_bp']==12 and len(r['path'])==1, r
cx.executemany('INSERT INTO edges VALUES (?,?,?,?)',[('u2','+','b1','+'),('u2','+','b2','+')])
r=g.walk(cx,('u0','+'),'RIGHT',100)
assert r['status']=='PARTIAL_GRAPH_CENSORED' and r['stop_reason']=='BRANCH', r
print('[PASS] F3B4A3 sqlite oriented unique-path/left-right/branch selftest')
PY
