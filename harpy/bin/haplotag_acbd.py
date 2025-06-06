#! /usr/bin/env python
"""Generates the BC_{ABCD}.txt files necessary to demultiplex Gen I haplotagging barcodes"""
import os
import sys
import argparse

parser = argparse.ArgumentParser(
    prog = 'haplotag_acbd.py',
    description ="Generates the BC_{ABCD}.txt files necessary to demultiplex Gen I haplotagging barcodes",
    usage = "haplotag_acbd.py output_directory",
    exit_on_error = False
    )
parser.add_argument("output_directory", type = str, help = "Directory to create barcode files")
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
outdir = args.output_directory.rstrip("/")
os.makedirs(outdir, exist_ok = True)

BX = {
    "A": ["ACGGAA", "CCAACA", "AGATCG", "TTCTCC", "TTCCTG", "TTCGGT", "TTGTGG", "TTGCCT", "TTGGTC", "TTACGC", "TTAGCG", "TCTTCG", "TCTCTC", "TCTGGA", "TCCACT", "TCGTAC", "TCGATG", "TCACAG", "TGTTGC", "TGTCCA", "TGTGTG", "TGCTAG", "TGCATC", "TGGAGT", "TGAGAC", "TATCGG", "TATGCC", "TACCAC", "TAGGAG", "CTTCGT", "CTTGCA", "CTCTGA", "CTCAAC", "CTGCTA", "CTGGAT", "CTAAGG", "CCTCAA", "CCTGTT", "CCATTC", "CGTTCT", "CGTAGA", "CGGTAA", "CGACTT", "CATACG", "CACTTG", "CACGAA", "CACAGT", "CAGATC", "CAACGA", "CAAGCT", "GTTCAC", "GTCGTA", "GTGTCA", "GTGAAG", "GTAACC", "GCTTGT", "GCCTAA", "GCACTA", "GCAGAT", "GGTGAA", "GGCAAT", "GGATGA", "GGAATG", "GATCCT", "GATAGC", "GACACA", "GAGCAA", "GAGGTT", "ATTCCG", "ATTGGC", "ATCGAG", "ACTACC", "ACCAGA", "ACGTCT", "ACACGT", "ACAGTG", "AGCTGT", "AGCCTA", "AGGTTC", "AGGCAT", "AGGACA", "AGAAGC", "AACGTC", "AAGCTG", "CGAGTA", "GAATCC", "GAATGG", "AAGTGC", "AAGAGG", "TACAGG", "CTGACT", "CTAGTC", "CCTAAG", "CCATAG", "CGTAAC", "CAATGC"],
    "C": ["GAAACG", "ACACCA", "TCGAGA", "TCCTTC", "CTGTTC", "GGTTTC", "TGGTTG", "CCTTTG", "GTCTTG", "CGCTTA", "GCGTTA", "TCGTCT", "CTCTCT", "GGATCT", "ACTTCC", "TACTCG", "ATGTCG", "CAGTCA", "TGCTGT", "CCATGT", "GTGTGT", "TAGTGC", "ATCTGC", "AGTTGG", "GACTGA", "CGGTAT", "GCCTAT", "CACTAC", "GAGTAG", "CGTCTT", "GCACTT", "TGACTC", "AACCTC", "CTACTG", "GATCTG", "AGGCTA", "CAACCT", "GTTCCT", "TTCCCA", "TCTCGT", "AGACGT", "TAACGG", "CTTCGA", "ACGCAT", "TTGCAC", "GAACAC", "AGTCAC", "ATCCAG", "CGACAA", "GCTCAA", "CACGTT", "GTAGTC", "TCAGTG", "AAGGTG", "ACCGTA", "TGTGCT", "TAAGCC", "CTAGCA", "GATGCA", "GAAGGT", "AATGGC", "TGAGGA", "ATGGGA", "CCTGAT", "AGCGAT", "ACAGAC", "CAAGAG", "GTTGAG", "CCGATT", "GGCATT", "GAGATC", "ACCACT", "AGAACC", "TCTACG", "CGTACA", "GTGACA", "TGTAGC", "CTAAGC", "TTCAGG", "CATAGG", "ACAAGG", "AGCAGA", "GTCAAC", "CTGAAG", "GTACGA", "TCCGAA", "TGGGAA", "TGCAAG", "AGGAAG", "AGGTAC", "ACTCTG", "GTCCTA", "AAGCCT", "TAGCCA", "AACCGT", "TGCCAA"],
    "B": ["AACGGA", "ACCAAC", "GAGATC", "CTTCTC", "GTTCCT", "TTTCGG", "GTTGTG", "TTTGCC", "CTTGGT", "CTTACG", "GTTAGC", "GTCTTC", "CTCTCT", "ATCTGG", "TTCCAC", "CTCGTA", "GTCGAT", "GTCACA", "CTGTTG", "ATGTCC", "GTGTGT", "GTGCTA", "CTGCAT", "TTGGAG", "CTGAGA", "GTATCG", "CTATGC", "CTACCA", "GTAGGA", "TCTTCG", "ACTTGC", "ACTCTG", "CCTCAA", "ACTGCT", "TCTGGA", "GCTAAG", "ACCTCA", "TCCTGT", "CCCATT", "TCGTTC", "ACGTAG", "ACGGTA", "TCGACT", "GCATAC", "GCACTT", "ACACGA", "TCACAG", "CCAGAT", "ACAACG", "TCAAGC", "CGTTCA", "AGTCGT", "AGTGTC", "GGTGAA", "CGTAAC", "TGCTTG", "AGCCTA", "AGCACT", "TGCAGA", "AGGTGA", "TGGCAA", "AGGATG", "GGGAAT", "TGATCC", "CGATAG", "AGACAC", "AGAGCA", "TGAGGT", "GATTCC", "CATTGG", "GATCGA", "CACTAC", "AACCAG", "TACGTC", "TACACG", "GACAGT", "TAGCTG", "AAGCCT", "CAGGTT", "TAGGCA", "AAGGAC", "CAGAAG", "CAACGT", "GAAGCT", "ACGAGT", "CGAATC", "GGAATG", "CAAGTG", "GAAGAG", "GTACAG", "TCTGAC", "CCTAGT", "GCCTAA", "GCCATA", "CCGTAA", "CCAATG"],
    "D": ["GGAAAC", "AACACC", "ATCGAG", "CTCCTT", "CCTGTT", "CGGTTT", "GTGGTT", "GCCTTT", "GGTCTT", "ACGCTT", "AGCGTT", "TTCGTC", "TCTCTC", "TGGATC", "CACTTC", "GTACTC", "GATGTC", "ACAGTC", "TTGCTG", "TCCATG", "TGTGTG", "CTAGTG", "CATCTG", "GAGTTG", "AGACTG", "TCGGTA", "TGCCTA", "CCACTA", "GGAGTA", "TCGTCT", "TGCACT", "CTGACT", "CAACCT", "GCTACT", "GGATCT", "AAGGCT", "TCAACC", "TGTTCC", "ATTCCC", "TTCTCG", "TAGACG", "GTAACG", "ACTTCG", "TACGCA", "CTTGCA", "CGAACA", "CAGTCA", "GATCCA", "ACGACA", "AGCTCA", "TCACGT", "CGTAGT", "GTCAGT", "GAAGGT", "AACCGT", "TTGTGC", "CTAAGC", "ACTAGC", "AGATGC", "TGAAGG", "CAATGG", "ATGAGG", "AATGGG", "TCCTGA", "TAGCGA", "CACAGA", "GCAAGA", "GGTTGA", "TCCGAT", "TGGCAT", "CGAGAT", "TACCAC", "CAGAAC", "GTCTAC", "ACGTAC", "AGTGAC", "CTGTAG", "CCTAAG", "GTTCAG", "GCATAG", "GACAAG", "AAGCAG", "CGTCAA", "GCTGAA", "AGTACG", "ATCCGA", "ATGGGA", "GTGCAA", "GAGGAA", "CAGGTA", "GACTCT", "AGTCCT", "TAAGCC", "ATAGCC", "TAACCG", "ATGCCA"]
}

for BC in ["A","C","B","D"]:
    with open(f"{outdir}/segment_{BC}.bc", "w", encoding="utf-8") as f:
        ID = [f"{BC}{number:02d}" for number in range(1, 97)]
        delim = ["	".join(tup) for tup in zip(ID, BX[BC])]
        _ = [f.write(f"{i}\n") for i in delim]