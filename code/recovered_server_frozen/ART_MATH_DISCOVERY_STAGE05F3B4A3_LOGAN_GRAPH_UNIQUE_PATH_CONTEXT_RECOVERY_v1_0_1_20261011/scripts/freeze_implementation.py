#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]); b4a2=Path(sys.argv[2]); b4a=Path(sys.argv[3])
def sha(p):
 h=hashlib.sha256(); h.update(Path(p).read_bytes()); return h.hexdigest()
g={
 'stage':'Stage05F3B4A3-implementation-freeze','gate':'PASS_F3B4A3_GRAPH_UNIQUE_PATH_RULE_FROZEN_PRE_GRAPH_RESULTS','gate_pass':True,
 'input_survivor_n':16,'input_boundary_anchor_n':31,'required_anchor_resolution':'UNIQUE','kmer_k':31,'unitig_overlap_bp':30,
 'graph_source':'Logan k=31 unitigs with BCALM2-style L:<+/->:<other id>:<+/-> links',
 'extension_rule':'for each missing boundary, orient the unique anchor unitig to the parent-contig boundary; walk only a single degree-1 oriented continuation; stop at branch/dead-end/cycle/missing link or when missing flank is recovered',
 'branch_choice_rule':'FORBIDDEN; no abundance, annotation, repeat, RT, ART, or downstream result may select among multiple graph branches',
 'target_flank_each_side_bp':10000,'repeat_analysis_performed_in_stage':False,'partner_analysis_performed_in_stage':False,'novelty_decision_performed':False,
 'b4a_gate_sha256':sha(b4a/'outputs/STAGE05F3B4A_GATE.json'),
 'b4a2_context_gate_sha256':sha(b4a2/'outputs/STAGE05F3B4A2_CONTEXT_AUDIT_GATE.json'),
 'b4a2_anchor_scan_sha256':sha(b4a2/'outputs/F3B4A2_LOGAN_UNITIG_BOUNDARY_ANCHOR_SCAN.tsv'),
 'implementation_sha256':sha(root/'scripts/graph_unique_path_recovery.py'),
 'claim_boundary':'Technical source-context reconstruction under a frozen unique-path graph rule only; no ART architecture or novelty claim.'
}
(root/'outputs').mkdir(exist_ok=True); (root/'outputs/STAGE05F3B4A3_IMPLEMENTATION_FREEZE.json').write_text(json.dumps(g,indent=2)+'\n')
print(json.dumps(g,indent=2))
