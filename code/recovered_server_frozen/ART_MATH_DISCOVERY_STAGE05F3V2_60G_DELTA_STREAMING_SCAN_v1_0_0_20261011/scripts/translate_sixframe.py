#!/usr/bin/env python3
import sys,argparse,re,json
ap=argparse.ArgumentParser(); ap.add_argument('--accession',required=True); ap.add_argument('--min-aa',type=int,default=45); ap.add_argument('--stats',required=True)
a=ap.parse_args()
# Genetic code 11 amino-acid mapping is identical to standard code for internal codons; no start-codon requirement is used.
code={
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
comp=str.maketrans('ACGTRYMKSWBDHVNacgtrymkswbdhvn','TGCAYRKMSWVHDBNtgcayrkmswvhdbn')
def rc(s): return s.translate(comp)[::-1].upper()
def safe(x): return re.sub(r'[^A-Za-z0-9_.:-]+','_',x)[:180]

def emit(contig,seq):
    global orfs,aa_total
    L=len(seq); cid=safe(contig)
    for strand,ss in [('+',seq.upper()),('-',rc(seq))]:
      for frame in (0,1,2):
        seg=[]; seg_codon0=0; codoni=0
        ncod=(len(ss)-frame)//3
        for j in range(ncod):
            cod=ss[frame+3*j:frame+3*j+3]
            aa=code.get(cod,'X')
            if aa=='*':
                if len(seg)>=a.min_aa:
                    r0=frame+3*seg_codon0; r1=frame+3*j
                    if strand=='+': s0,e0=r0,r1
                    else: s0,e0=L-r1,L-r0
                    q=f'{a.accession}|{cid}|{strand}|f{frame}|{s0}-{e0}|aa{len(seg)}|o{orfs}'
                    print('>'+q); print(''.join(seg)); orfs+=1; aa_total+=len(seg)
                seg=[]; seg_codon0=j+1
            else: seg.append(aa)
        if len(seg)>=a.min_aa:
            r0=frame+3*seg_codon0; r1=frame+3*ncod
            if strand=='+': s0,e0=r0,r1
            else: s0,e0=L-r1,L-r0
            q=f'{a.accession}|{cid}|{strand}|f{frame}|{s0}-{e0}|aa{len(seg)}|o{orfs}'
            print('>'+q); print(''.join(seg)); orfs+=1; aa_total+=len(seg)

name=None; chunks=[]; contigs=0; nt_total=0; orfs=0; aa_total=0
for raw in sys.stdin:
    line=raw.strip()
    if not line: continue
    if line.startswith('>'):
        if name is not None:
            s=''.join(chunks); contigs+=1; nt_total+=len(s); emit(name,s)
        name=line[1:].split()[0]; chunks=[]
    else: chunks.append(line)
if name is not None:
    s=''.join(chunks); contigs+=1; nt_total+=len(s); emit(name,s)
open(a.stats,'w').write(json.dumps({'accession':a.accession,'contigs':contigs,'nt_bp':nt_total,'orf_n':orfs,'aa_total':aa_total},indent=2)+'\n')
