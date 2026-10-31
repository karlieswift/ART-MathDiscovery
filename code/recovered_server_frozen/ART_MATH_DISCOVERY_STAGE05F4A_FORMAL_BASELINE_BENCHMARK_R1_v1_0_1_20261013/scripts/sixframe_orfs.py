#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
CODON={
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
COMP=str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(COMP)[::-1]
def translate(s): return ''.join(CODON.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
def sixframe(seq,min_aa=45):
    n=len(seq)
    for strand,dna in [('+',seq),('-',rc(seq))]:
        for frame in range(3):
            aa=translate(dna[frame:])
            start=0
            for j in range(len(aa)+1):
                if j==len(aa) or aa[j]=='*':
                    pep=aa[start:j]
                    if len(pep)>=min_aa:
                        nt0=frame+start*3; nt1=frame+j*3
                        if strand=='+': g0,g1=nt0,nt1
                        else: g0,g1=n-nt1,n-nt0
                        yield strand,frame,start,j,g0,g1,pep
                    start=j+1
