#! /usr/bin/env python
"""Generates a text file listing the haplotagging ACBD barcodes"""
import sys
import argparse
from itertools import product

parser = argparse.ArgumentParser(
    prog = 'haplotag_barcodes.py',
    description ="Prints haplotagging linked-read barcodes to stdout",
    usage = "haplotag_barcodes.py [-n] > barcodes.txt"
    )
parser.add_argument("-n", "--number", default=96**4, type = int, help = "How many barcodes to print (min: 1, max: 84934656)")
args = parser.parse_args()

if args.number < 1 or args.number > 96**4:
    parser.error("--number must be between 1 and 96^4 (84934656)")
# subtract 1 b/c of python 0-indexing
args.number -= 1

BX = {
    "A": ["ACGGAA", "CCAACA", "AGATCG", "TTCTCC", "TTCCTG", "TTCGGT", "TTGTGG", "TTGCCT", "TTGGTC", "TTACGC", "TTAGCG", "TCTTCG", "TCTCTC", "TCTGGA", "TCCACT", "TCGTAC", "TCGATG", "TCACAG", "TGTTGC", "TGTCCA", "TGTGTG", "TGCTAG", "TGCATC", "TGGAGT", "TGAGAC", "TATCGG", "TATGCC", "TACCAC", "TAGGAG", "CTTCGT", "CTTGCA", "CTCTGA", "CTCAAC", "CTGCTA", "CTGGAT", "CTAAGG", "CCTCAA", "CCTGTT", "CCATTC", "CGTTCT", "CGTAGA", "CGGTAA", "CGACTT", "CATACG", "CACTTG", "CACGAA", "CACAGT", "CAGATC", "CAACGA", "CAAGCT", "GTTCAC", "GTCGTA", "GTGTCA", "GTGAAG", "GTAACC", "GCTTGT", "GCCTAA", "GCACTA", "GCAGAT", "GGTGAA", "GGCAAT", "GGATGA", "GGAATG", "GATCCT", "GATAGC", "GACACA", "GAGCAA", "GAGGTT", "ATTCCG", "ATTGGC", "ATCGAG", "ACTACC", "ACCAGA", "ACGTCT", "ACACGT", "ACAGTG", "AGCTGT", "AGCCTA", "AGGTTC", "AGGCAT", "AGGACA", "AGAAGC", "AACGTC", "AAGCTG", "CGAGTA", "GAATCC", "GAATGG", "AAGTGC", "AAGAGG", "TACAGG", "CTGACT", "CTAGTC", "CCTAAG", "CCATAG", "CGTAAC", "CAATGC"],
    "C": ["GAAACG", "ACACCA", "TCGAGA", "TCCTTC", "CTGTTC", "GGTTTC", "TGGTTG", "CCTTTG", "GTCTTG", "CGCTTA", "GCGTTA", "TCGTCT", "CTCTCT", "GGATCT", "ACTTCC", "TACTCG", "ATGTCG", "CAGTCA", "TGCTGT", "CCATGT", "GTGTGT", "TAGTGC", "ATCTGC", "AGTTGG", "GACTGA", "CGGTAT", "GCCTAT", "CACTAC", "GAGTAG", "CGTCTT", "GCACTT", "TGACTC", "AACCTC", "CTACTG", "GATCTG", "AGGCTA", "CAACCT", "GTTCCT", "TTCCCA", "TCTCGT", "AGACGT", "TAACGG", "CTTCGA", "ACGCAT", "TTGCAC", "GAACAC", "AGTCAC", "ATCCAG", "CGACAA", "GCTCAA", "CACGTT", "GTAGTC", "TCAGTG", "AAGGTG", "ACCGTA", "TGTGCT", "TAAGCC", "CTAGCA", "GATGCA", "GAAGGT", "AATGGC", "TGAGGA", "ATGGGA", "CCTGAT", "AGCGAT", "ACAGAC", "CAAGAG", "GTTGAG", "CCGATT", "GGCATT", "GAGATC", "ACCACT", "AGAACC", "TCTACG", "CGTACA", "GTGACA", "TGTAGC", "CTAAGC", "TTCAGG", "CATAGG", "ACAAGG", "AGCAGA", "GTCAAC", "CTGAAG", "GTACGA", "TCCGAA", "TGGGAA", "TGCAAG", "AGGAAG", "AGGTAC", "ACTCTG", "GTCCTA", "AAGCCT", "TAGCCA", "AACCGT", "TGCCAA"],
    "B": ["AACGGA", "ACCAAC", "GAGATC", "CTTCTC", "GTTCCT", "TTTCGG", "GTTGTG", "TTTGCC", "CTTGGT", "CTTACG", "GTTAGC", "GTCTTC", "CTCTCT", "ATCTGG", "TTCCAC", "CTCGTA", "GTCGAT", "GTCACA", "CTGTTG", "ATGTCC", "GTGTGT", "GTGCTA", "CTGCAT", "TTGGAG", "CTGAGA", "GTATCG", "CTATGC", "CTACCA", "GTAGGA", "TCTTCG", "ACTTGC", "ACTCTG", "CCTCAA", "ACTGCT", "TCTGGA", "GCTAAG", "ACCTCA", "TCCTGT", "CCCATT", "TCGTTC", "ACGTAG", "ACGGTA", "TCGACT", "GCATAC", "GCACTT", "ACACGA", "TCACAG", "CCAGAT", "ACAACG", "TCAAGC", "CGTTCA", "AGTCGT", "AGTGTC", "GGTGAA", "CGTAAC", "TGCTTG", "AGCCTA", "AGCACT", "TGCAGA", "AGGTGA", "TGGCAA", "AGGATG", "GGGAAT", "TGATCC", "CGATAG", "AGACAC", "AGAGCA", "TGAGGT", "GATTCC", "CATTGG", "GATCGA", "CACTAC", "AACCAG", "TACGTC", "TACACG", "GACAGT", "TAGCTG", "AAGCCT", "CAGGTT", "TAGGCA", "AAGGAC", "CAGAAG", "CAACGT", "GAAGCT", "ACGAGT", "CGAATC", "GGAATG", "CAAGTG", "GAAGAG", "GTACAG", "TCTGAC", "CCTAGT", "GCCTAA", "GCCATA", "CCGTAA", "CCAATG"],
    "D": ["GGAAAC", "AACACC", "ATCGAG", "CTCCTT", "CCTGTT", "CGGTTT", "GTGGTT", "GCCTTT", "GGTCTT", "ACGCTT", "AGCGTT", "TTCGTC", "TCTCTC", "TGGATC", "CACTTC", "GTACTC", "GATGTC", "ACAGTC", "TTGCTG", "TCCATG", "TGTGTG", "CTAGTG", "CATCTG", "GAGTTG", "AGACTG", "TCGGTA", "TGCCTA", "CCACTA", "GGAGTA", "TCGTCT", "TGCACT", "CTGACT", "CAACCT", "GCTACT", "GGATCT", "AAGGCT", "TCAACC", "TGTTCC", "ATTCCC", "TTCTCG", "TAGACG", "GTAACG", "ACTTCG", "TACGCA", "CTTGCA", "CGAACA", "CAGTCA", "GATCCA", "ACGACA", "AGCTCA", "TCACGT", "CGTAGT", "GTCAGT", "GAAGGT", "AACCGT", "TTGTGC", "CTAAGC", "ACTAGC", "AGATGC", "TGAAGG", "CAATGG", "ATGAGG", "AATGGG", "TCCTGA", "TAGCGA", "CACAGA", "GCAAGA", "GGTTGA", "TCCGAT", "TGGCAT", "CGAGAT", "TACCAC", "CAGAAC", "GTCTAC", "ACGTAC", "AGTGAC", "CTGTAG", "CCTAAG", "GTTCAG", "GCATAG", "GACAAG", "AAGCAG", "CGTCAA", "GCTGAA", "AGTACG", "ATCCGA", "ATGGGA", "GTGCAA", "GAGGAA", "CAGGTA", "GACTCT", "AGTCCT", "TAAGCC", "ATAGCC", "TAACCG", "ATGCCA"]
}
   
bc_generator = product(BX["A"], BX["C"], BX["B"], BX["D"])
for i,BC in enumerate(bc_generator):
    sys.stdout.write("".join(BC) + "\n")
    if i >= args.number:
        break