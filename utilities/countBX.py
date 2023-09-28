#! /usr/bin/env python3

import pysam
import re
import argparse
import sys

parser = argparse.ArgumentParser(
    prog = 'countBX.py',
    description = 'Count the number of valid Haplotag BX tags in a FASTQ file.',
    usage = "countBX.py fastqfile",
    exit_on_error = False
    )

parser.add_argument("fastqfile", help = "Input FASTQ file.")

# parser.add_argument("outfile", help = "Output bam file. A corresponding index file will be created for it.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

n_reads = 0
n_bx = 0
n_valid = 0
haplotag = re.compile('A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}')
invalid = re.compile('[AaBbCcDd]00')
inv_dict = {
    "A" : 0,
    "B" : 0,
    "C" : 0,
    "D" : 0
}
with pysam.FastxFile(args.fastqfile) as fh:
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
                        inv_dict[i[0]] += 1
                    continue
                n_valid += 1
            

print(f"totalReads\t{n_reads}")
print(f"bxTagCount\t{n_bx}")
print(f"bxValid\t{n_valid}")
print(f"bxInvalid\t{n_bx - n_valid}")
print("A00\t",str(inv_dict["A"]))
print("C00\t",str(inv_dict["C"]))
print("B00\t",str(inv_dict["B"]))
print("D00\t",str(inv_dict["D"]))