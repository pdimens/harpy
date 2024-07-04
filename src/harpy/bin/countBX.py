#! /usr/bin/env python

import re
import sys
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'countBX.py',
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

n_reads = 0
n_bx = 0
n_valid = 0
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
        n_reads += 1
        comments = entry.comment.split()
        # looking for a comment that starts as whitespace + BX:Z:
        bxtag_idx = [i for i,j in enumerate(comments) if j.startswith("BX:Z:")]
        #if 'BX:Z:' in entry.comment:
        if bxtag_idx:
            n_bx += 1
            beadtag_full = comments[bxtag_idx[0]]
            beadtag = beadtag_full[5:]
            if bool(haplotag.match(beadtag)):
                inv = re.findall(invalid, beadtag)
                if inv:
                    for i in inv:
                    #    if i[0] in inv_dict:
                        inv_dict[i[0]] += 1
                    #   else:
                    #   inv_dict[i[0]] = 1
                    continue
                n_valid += 1

print(f"totalReads\t{n_reads}", file = sys.stdout)
print(f"bxTagCount\t{n_bx}", file = sys.stdout)
print(f"bxValid\t{n_valid}", file = sys.stdout)
print(f"bxInvalid\t{n_bx - n_valid}", file = sys.stdout)
print("A00\t",str(inv_dict["A"]), file = sys.stdout)
print("C00\t",str(inv_dict["C"]), file = sys.stdout)
print("B00\t",str(inv_dict["B"]), file = sys.stdout)
print("D00\t",str(inv_dict["D"]), file = sys.stdout)
