from pathlib import Path
import tempfile,sys,unittest
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from scan_groupwise import metadata,starts_for,BLOCKS
from io_utils import read_json
import reference_math as ref

class ReplayTests(unittest.TestCase):
    def test_complete_strata_reproduce_full_matrix_selection(self):
        df=metadata(2200,'synthetic',[128,256,512,1024,2048],0.25)
        rng=np.random.default_rng(602905);X=rng.integers(0,6,size=(len(df),5468)).astype('float32')
        X[:,0]=1;X[:,1]=0;X[0,1]=1e6;X[:,4]*=1e-14
        for name,cols in BLOCKS.items():
            full=np.zeros(len(df));grouped=np.zeros(len(df))
            for (L,strand),ii0 in df.groupby(['window_len','strand'],sort=True).groups.items():
                ii=np.asarray(list(ii0),int)
                a,k=ref.robust_score(X[np.ix_(ii,cols)])
                # Simulate the bounded extraction file without changing the
                # complete stratum or the original advanced indexing layout.
                with tempfile.TemporaryDirectory() as td:
                    p=Path(td)/'group.npy';M=np.lib.format.open_memmap(p,mode='w+',shape=(len(ii),5468),dtype='float32')
                    for base in range(0,len(ii),7):M[base:base+7]=X[ii[base:base+7]]
                    b,kb=ref.robust_score(M[np.ix_(np.arange(len(ii)),cols)])
                    self.assertEqual(k,kb);np.testing.assert_array_equal(a,b)
                    full[ii]=a;grouped[ii]=b
                    del M
            pd.testing.assert_frame_equal(ref.collapse_oriented_q(df,full,name),ref.collapse_oriented_q(df,grouped,name),check_exact=True)
    def test_chunk_calibration_is_not_equivalent(self):
        X=np.array([[0],[0],[0],[0],[10],[10],[10],[10]],dtype=np.float32)
        full,_=ref.robust_score(X)
        chunks=np.r_[ref.robust_score(X[:4])[0],ref.robust_score(X[4:])[0]]
        self.assertFalse(np.array_equal(full,chunks))
    def test_ties_keep_frozen_plus_precedence(self):
        df=metadata(256,'s',[128],0.25)
        b=ref.collapse_oriented_q(df,np.ones(len(df)),'a')
        self.assertTrue((b.a_best_strand=='+').all())
    def test_zero_scale_and_nonfinite(self):
        x,k=ref.robust_score(np.ones((20,11),dtype='float32'))
        self.assertEqual(k,0);self.assertTrue((x==0).all())
        with self.assertRaises(ValueError):ref.robust_score([[float('nan')]])
    def test_tail_and_short_contigs(self):
        self.assertEqual(starts_for(127,128,.25),[])
        self.assertEqual(starts_for(128,128,.25),[0])
        self.assertEqual(starts_for(161,128,.25),[0,32,33])
    def test_strict_overlap_does_not_merge_touching(self):
        x=pd.DataFrame({'start0':[0,10,11],'end0':[10,12,20]})
        self.assertEqual(ref.merge_overlapping(x),[[0,10],[10,20]])
    def test_frozen_15_contig_resource_arithmetic(self):
        p=read_json(Path(__file__).resolve().parents[1]/'config/POLICY.json')
        counts=[];peaks=[]
        for r in p['contigs']:
            ss=[len(starts_for(r['length_bp'],L,.25)) for L in [128,256,512,1024,2048,4096]]
            counts.append(2*sum(ss));peaks.append(max(ss)*5468*4)
        self.assertEqual(sum(counts),2542244)
        self.assertEqual(sum(r['length_bp'] for r in p['contigs']),20664383)
        self.assertLess(max(peaks),4_000_000_000)
        self.assertEqual(len(p['contigs']),15)
        self.assertEqual(len({r['sequence_sha256'] for r in p['contigs']}),14)

    def test_synthetic_group_worker_integration(self):
        import hashlib,json
        import scan_groupwise as scan
        from io_utils import single_fasta,digest
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);eng=root/'engine'
            for d in ('scripts','config','frozen_sources'):(eng/d).mkdir(parents=True)
            mathcode=Path(__file__).with_name('reference_math.py').read_text()
            mathcode += "\nfrom io_utils import read_json, single_fasta\ndef parse_fasta(p):return single_fasta(p)\ndef revcomp(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]\ndef recurrence_scores(seq,windows,*args):\n n=len(windows);a=(windows.start0.to_numpy()%13).astype(float);return a,np.full(n,4,dtype=np.int16),np.full(n,8,dtype=np.int16),np.full(n,2,dtype=np.int32)\n"
            (eng/'scripts/common.py').write_text(mathcode)
            fixture = """from common import *
import numpy as np,hashlib
def init_acnv(src):pass
def init_irr(src,vocab):pass
def feat(seq,width):
 seed=int.from_bytes(hashlib.sha256(seq.encode()).digest()[:4],'little')
 return np.random.default_rng(seed).integers(0,7,size=width).astype('float32')
def acnv_one(seq):return feat(seq,5468)
def irr_one(item):return item[0],feat(item[1],3570)
def merge_regions(selected,d):
 fams=[('ACNV_FAMILY','q_acnv_family'),('IRREGULARITY_V3_D20','q_irr'),('GENERIC_RECURRENCE','q_recurrence')]
 rows=[]
 for start,end in merge_overlapping(selected):
  r={'region_rank':len(rows)+1,'region_id':'R%04d'%(len(rows)+1),'start0':start,'end0':end,'span_bp':end-start}
  support=[]
  for fam,col in fams:
   v=selected[(selected.start0<end)&(selected.end0>start)&(selected.family==fam)]
   if len(v):support.append(fam)
  r['support_families']=';'.join(support);rows.append(r)
 return pd.DataFrame(rows)
"""
            (eng/'scripts/02_primary_scan.py').write_text(fixture)
            pol={'window_lengths_bp':[128,256], 'stride_fraction':0.25,'proposal_q_threshold':0.2,
                 'recurrence':{'k_values':[4,6,8],'min_lag_bp':8,'max_lag_bp':512,'min_valid_pairs':32,'min_recurrence_match_count':2}}
            (eng/'config/ENGINE_POLICY.json').write_text(json.dumps(pol))
            fasta=root/'test.fa';seq=''.join(np.random.default_rng(50).choice(list('ACGT'),1024));fasta.write_text('>SYNTHETIC\n'+seq+'\n')
            import contextlib,io
            with contextlib.redirect_stdout(io.StringIO()):g=scan.run(eng,fasta,'SYNTHETIC',root/'out',root/'work',2)
            cm,pr,fams=scan.load_engine(eng);df=scan.metadata(len(seq),'SYNTHETIC',pol['window_lengths_bp'],pol['stride_fraction'])
            zz=[seq[r.start0:r.end0] if r.strand=='+' else cm.revcomp(seq[r.start0:r.end0]) for r in df.itertuples()]
            ac=np.stack([pr.acnv_one(z) for z in zz]);aux=np.stack([pr.irr_one((i,z))[1] for i,z in enumerate(zz)])
            scores={name:np.zeros(len(df)) for name in list(scan.BLOCKS)+['irr']}
            for _,ii0 in df.groupby(['window_len','strand'],sort=True).groups.items():
                ii=np.asarray(list(ii0),int)
                for name,col in scan.BLOCKS.items():scores[name][ii]=cm.robust_score(ac[np.ix_(ii,col)])[0]
                scores['irr'][ii]=cm.robust_score(aux[ii])[0]
            tables=scan.assemble_boards(cm,pr,fams,seq,df,scores,pol)
            for name,tab in zip(scan.BOARDS,tables):
                import gzip
                actual=(root/'out'/name).read_bytes()
                if name.endswith('.gz'):actual=gzip.decompress(actual)
                self.assertEqual(tab.to_csv(sep='\t',index=False).encode(),actual)
            with contextlib.redirect_stdout(io.StringIO()):g2=scan.run(eng,fasta,'SYNTHETIC',root/'out',root/'work',2)
            self.assertEqual(g['boards'],g2['boards'])
            self.assertFalse(list((root/'work').glob('*.npy')))

if __name__=='__main__':unittest.main()
