#! /usr/bin/env python3

"""Creates a FASTA file with 2 contigs of 200k random nucleotides""" 
from itertools import product
from random import choices

nuc = "ATCG"
with open("test.hap1.fa", "w") as fa:
    for i in [1,2]:
        _ = fa.write(f">Contig{i}\n" + "".join(choices(nuc, k = 200000)) + "\n")

with open("test.hap2.fa", "w") as fa:
    for i in [1,2]:
        _ = fa.write(f">Contig{i}\n" + "".join(choices(nuc, k = 200000)) + "\n")

with open("nucleotides.bc", "w") as bc:
    bc_gen = product(*["ATCG" for i in range(18)])
    for i in range(1000000):
        _bc = "".join(next(bc_gen))
        bc.write(_bc + "\n")

with open("combinatorial.bc", "w") as bc:
    bc_gen = product(*["ATCG" for i in range(6)])
    for i in range(96):
        _bc = "".join(next(bc_gen))
        bc.write(_bc + "\n")