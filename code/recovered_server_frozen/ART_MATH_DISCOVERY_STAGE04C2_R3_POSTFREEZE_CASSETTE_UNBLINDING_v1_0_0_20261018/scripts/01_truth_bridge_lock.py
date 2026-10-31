#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
from common import sha256_file

EXPECTED_BRIDGE='8f3e1cc6d99927159a482d250cd8db3d09c522f1d1ae89c383a5cf7a7c99c780'
EXPECTED_GATE='1fb15294b73a9a5bbf8769abf6dac5745408adf4520bc7a86b3415870348bab3'
EXPECTED_INPUT_LOCK='2d60695a3ab7ce0b0b09ff3a68236201789918f78ff29c035f3bbe56fa867f21'

def fail(x): raise SystemExit('HOLD_STAGE04C2_R3_TRUTH_BRIDGE_LOCK: '+x)

ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--r2c',required=True); ap.add_argument('--ug',required=True); a=ap.parse_args()
root=Path(a.root).resolve(); r2c=Path(a.r2c).resolve(); ug=Path(a.ug).resolve(); out=root/'outputs'
pre=json.loads((out/'STAGE04C2_R3_SCORE_ARCHIVE_PREFLIGHT_GATE.json').read_text())
if not pre.get('gate_pass'): fail('score archive preflight did not pass')

bridge=ug/'outputs/stage04c2b_r3/UG27_FIVE_SYSTEM_CONSTRUCT_NATIVE_BRIDGE.tsv'
gatep=ug/'outputs/stage04c2b_r3/GATE_STAGE04C2B_R3_CONSTRUCT_NATIVE_BRIDGE.json'
lockp=r2c/'outputs/CONTIG_INPUT_LOCK.tsv'
if not bridge.is_file() or sha256_file(bridge)!=EXPECTED_BRIDGE: fail('bridge hash mismatch')
if not gatep.is_file() or sha256_file(gatep)!=EXPECTED_GATE: fail('bridge gate hash mismatch')
if not lockp.is_file() or sha256_file(lockp)!=EXPECTED_INPUT_LOCK: fail('contig input lock hash mismatch')

G=json.loads(gatep.read_text())
checks={
 'gate':G.get('gate')=='PASS_STAGE04C2B_R3_UG27_FIVE_CONSTRUCTS_BRIDGED_TO_COGNATE_NATIVE_LOCI',
 'gate_pass':G.get('gate_pass') is True,
 'expected5':int(G.get('expected_systems_n',-1))==5,
 'bridges5':int(G.get('constructs_with_native_bridge_n',-1))==5,
 'single5':int(G.get('single_cognate_native_hsp_n',-1))==5,
 'complete':G.get('construct_to_native_locus_bridge_complete') is True,
 'retuning_false':G.get('retuning_authorized') is False,
 'score_recomputed_false':G.get('score_recomputed') is False,
}
if not all(checks.values()): fail('bridge semantic checks failed: '+repr(checks))

B=list(csv.DictReader(bridge.open(),delimiter='\t'))
L=list(csv.DictReader(lockp.open(),delimiter='\t'))
if len(B)!=5: fail(f'expected 5 bridge rows, got {len(B)}')
by_src={}
for x in L: by_src.setdefault(x['source_contig_id'],[]).append(x)
rows=[]
for b in B:
    hits=by_src.get(b['sseqid'],[])
    if len(hits)!=1: fail(f'bridge subject does not map uniquely to frozen scoring contig: {b["sseqid"]} -> {len(hits)}')
    x=hits[0]
    s1=int(b['sstart']); s2=int(b['send']); start0=min(s1,s2)-1; end0=max(s1,s2)
    span=end0-start0
    if span!=int(b['aln_len']): fail(f'native HSP span != aln_len for {b["system"]}: {span} vs {b["aln_len"]}')
    if int(b['slen'])!=int(x['length_bp']): fail(f'slen != frozen contig length for {b["system"]}')
    if start0<0 or end0>int(x['length_bp']): fail(f'native interval out of bounds for {b["system"]}')
    orient='+' if s2>=s1 else '-'
    if orient!=b['orientation']: fail(f'orientation mismatch for {b["system"]}')
    rows.append({
      'system':b['system'],'assembly_id':x['assembly_id'],'scoring_id':x['scoring_id'],
      'source_contig_id':x['source_contig_id'],'contig_length_bp':x['length_bp'],
      'qseqid':b['qseqid'],'pident':b['pident'],'aln_len':b['aln_len'],'mismatch':b['mismatch'],'gapopen':b['gapopen'],
      'orientation':b['orientation'],'source_sstart_1based':b['sstart'],'source_send_1based':b['send'],
      'cassette_start0':start0,'cassette_end0':end0,'cassette_span_bp':span,
      'qstart':b['qstart'],'qend':b['qend'],'qlen':b['qlen'],'construct_native_fraction':b['construct_native_fraction'],
    })
rows.sort(key=lambda z:(z['assembly_id'],z['system']))
fields=list(rows[0])
outp=out/'UG27_R3_CASSETTE_TRUTH_LOCK.tsv'
with outp.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
D={
 'stage':'Stage04C2-R3-truth-lock',
 'gate':'PASS_STAGE04C2_R3_FIVE_FROZEN_CONSTRUCT_NATIVE_CASSETTES_MAPPED_TO_FROZEN_SCORING_CONTIGS',
 'gate_pass':True,'n_systems':5,'n_unique_scoring_contigs':len({r['scoring_id'] for r in rows}),
 'bridge_sha256':EXPECTED_BRIDGE,'bridge_gate_sha256':EXPECTED_GATE,'contig_input_lock_sha256':EXPECTED_INPUT_LOCK,
 'truth_type':'WT construct-to-native long-HSP cassette interval',
 'coordinate_system':'converted from BLAST 1-based inclusive subject coordinates to 0-based half-open',
 'exact_ncrna_array_truth':False,
 'truth_opened_after_r2d_archive_lock':True,
 'retuning_authorized':False,'score_recomputed':False,'q_recomputed':False,'science_change':False,
 'truth_lock_sha256':sha256_file(outp),
}
(out/'STAGE04C2_R3_TRUTH_LOCK_GATE.json').write_text(json.dumps(D,indent=2,sort_keys=True)+'\n')
print(outp.read_text())
print(json.dumps(D,indent=2,sort_keys=True))
