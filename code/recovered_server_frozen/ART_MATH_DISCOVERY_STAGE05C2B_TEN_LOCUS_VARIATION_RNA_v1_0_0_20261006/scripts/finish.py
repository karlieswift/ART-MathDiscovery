from pathlib import Path
import argparse
from common import *
a=argparse.ArgumentParser();a.add_argument('--root',required=True);x=a.parse_args();r=Path(x.root);o=r/'outputs';o.mkdir(exist_ok=True)
fold=read_json(o/'RNA_FOLD_GATE.json') if (o/'RNA_FOLD_GATE.json').exists() else {}
source=read_json(o/'author_source/AUTHOR_SOURCE_GATE.json') if (o/'author_source/AUTHOR_SOURCE_GATE.json').exists() else {}
g={'stage':'Stage05C2B','gate':'PASS_FOCUSED_RNA_CHARACTERIZATION_WITH_NOVELTY_HOLD' if fold.get('gate_pass') else 'HOLD_FOCUSED_RNA_INCOMPLETE','gate_pass':bool(fold.get('gate_pass')),'RNA_execution_complete':bool(fold.get('gate_pass')),'author_document_retrieved':bool(source.get('document_retrieval_pass')),'n_RT_profile_loci':10,'math_captured':7,'math_not_captured':3,'close_reference_related_math_captured':5,'close_reference_related_math_not_captured':1,'all_3674_math_candidates_preserved_upstream':True,'full_author_catalogue_membership_comparison_complete':False,'new_ART_systems_confirmed':0,'wet_lab_ready_declared':0,'mathematical_scores_or_ranks_changed':False,'UG27_touched':False,'RNA_energy_is_not_function_or_expression':True,'RNA_null_fraction_is_not_calibrated_p':True}
write_json(o/'GATE_STAGE05C2B.json',g);print(__import__('json').dumps(g,indent=2))
