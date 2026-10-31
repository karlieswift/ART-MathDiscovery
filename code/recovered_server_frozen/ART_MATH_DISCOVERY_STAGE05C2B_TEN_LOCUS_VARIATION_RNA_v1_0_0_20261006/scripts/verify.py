from pathlib import Path
import argparse
from common import validate_manifest
a=argparse.ArgumentParser();a.add_argument('--root',required=True);x=a.parse_args();n=validate_manifest(Path(x.root));print('[PASS] Package files validated:',n)
