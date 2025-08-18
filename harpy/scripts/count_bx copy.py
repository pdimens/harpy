#! /usr/bin/env python
"""parse a fastq file to count BX stats"""
import os
import re
import sys
import argparse
import pysam

def main():
    parser = argparse.ArgumentParser(
        prog = 'count_bx',
        description =
        """
        Parses a FASTQ file to count: total sequences, total number of BX tags,
        number of valid BX tags, number of invalid BX tags, number of
        invalid BX tag segments (i.e. A00, C00, B00, D00)
        """,
        usage = "count_bx input.fastq > output.txt",
        exit_on_error = False
        )

    parser.add_argument('input', help = "Input fastq file. Can be gzipped.")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not os.path.exists(args.input):
        parser.error(f"{args.input} was not found")

    N_READS = 0
    N_BX = 0
    N_VALID = 0
    haplotag = re.compile('A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}')
    invalid = re.compile('[AaBbCcDd]00')
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

    sys.stdout.write(f"totalReads\t{N_READS}\n")
    sys.stdout.write(f"bxTagCount\t{N_BX}\n")
    sys.stdout.write(f"bxValid\t{N_VALID}\n")
    sys.stdout.write(f"bxInvalid\t{N_BX - N_VALID}\n")
    sys.stdout.write(f"A00\t{inv_dict['A']}\n")
    sys.stdout.write(f"C00\t{inv_dict['C']}\n")
    sys.stdout.write(f"B00\t{inv_dict['B']}\n")
    sys.stdout.write(f"D00\t{inv_dict['D']}\n")
