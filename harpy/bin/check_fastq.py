#! /usr/bin/env python
"""parse a fastq file to do BX checks"""
import re
import os
import sys
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'check_fastq.py',
    description =
    """
    Parses a FASTQ file to check if any sequences don't conform to the SAM spec,
    whether BX:Z: is the last tag in the record, and the counts of: total reads,
    reads without BX:Z: tag, reads with incorrect BX:Z: tag.
    """,
    usage = "check_bam.py input.bam > output.txt",
    exit_on_error = False
    )

parser.add_argument('input', help = "Input fastq file. Can be gzipped.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
if not os.path.exists(args.input):
    parser.error(f"{args.input} was not found")

fq_in = args.input

#bxz = re.compile('BX:Z:')
samspec = re.compile('[A-Z][A-Z]:[AifZHB]:')
haplotag = re.compile('A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]')
bxlast = re.compile('BX:Z:A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]$')

with pysam.FastxFile(fq_in, persist=False) as fh:
    N_READS    = 0
    NO_BX       = 0
    BAD_BX      = 0
    BAD_SAM_SPEC = 0
    BX_NOT_LAST  = 0
    for entry in fh:
        N_READS += 1
        # look for BX:Z: tag
        if 'BX:Z:' in entry.comment:
            # if AXXCXXBXXDXX format isnt found, it's a bad BX
            if not haplotag.search(entry.comment):
                BAD_BX += 1
            splithead = entry.comment.split()
            for i in splithead:
                # if comments dont start with TAG:TYPE:, invalid SAM spec
                if not samspec.match(i):
                    BAD_SAM_SPEC += 1
            # if the BX:Z: isn't at the end, add to BX_NOT_LAST
            if splithead[-1].startswith('BX:Z'):
                pass
            else:
                BX_NOT_LAST += 1
        else:
            # missing BX:Z: tag
            NO_BX += 1

values = [str(i) for i in [os.path.basename(fq_in), N_READS, NO_BX, BAD_BX, BAD_SAM_SPEC, BX_NOT_LAST]]
sys.stdout.write("\t".join(values) + "\n")
