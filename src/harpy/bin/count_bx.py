#! /usr/bin/env python
"""parse a fastq file to count BX stats"""
import re
import sys
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'count_bx.py',
    description =
    """
    Parses a FASTQ file to count: total sequences, total number of BX tags,
    number of valid haplotagging BX tags, number of invalid BX tags, number of
    invalid BX tag segments (i.e. A00, C00, B00, D00)
    """,
    usage = "countBX.py input.fastq > output.txt",
    exit_on_error = False
    )

parser.add_argument('input', help = "Input fastq file. Can be gzipped.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

N_READS = 0
N_BX = 0
N_VALID = 0
# haplotag = re.compile("([A-Z]\d{2,}){3,}")
haplotag = re.compile('A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}')
# invalid = re.compile('[A-Z]00')
invalid = re.compile('[AaBbCcDd]00')
# inv_dict = {}
inv_dict = {
    "A" : 0,
    "B" : 0,
    "C" : 0,
    "D" : 0
}
with pysam.FastxFile(args.input) as fh:
    for entry in fh:
        N_READS += 1
        comments = entry.comment.split()
        # looking for a comment that starts as whitespace + BX:Z:
        bxtag_idx = [i for i,j in enumerate(comments) if j.startswith("BX:Z:")]
        #if 'BX:Z:' in entry.comment:
        if bxtag_idx:
            N_BX += 1
            beadtag_full = comments[bxtag_idx[0]]
            beadtag = beadtag_full[5:]
            if bool(haplotag.match(beadtag)):
                inv = re.findall(invalid, beadtag)
                if inv:
                    for i in inv:
                        inv_dict[i[0]] += 1
                    continue
                N_VALID += 1

print(f"totalReads\t{N_READS}", file = sys.stdout)
print(f"bxTagCount\t{N_BX}", file = sys.stdout)
print(f"bxValid\t{N_VALID}", file = sys.stdout)
print(f"bxInvalid\t{N_BX - N_VALID}", file = sys.stdout)
print("A00\t",str(inv_dict["A"]), file = sys.stdout)
print("C00\t",str(inv_dict["C"]), file = sys.stdout)
print("B00\t",str(inv_dict["B"]), file = sys.stdout)
print("D00\t",str(inv_dict["D"]), file = sys.stdout)
