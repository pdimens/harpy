#! /usr/bin/env python

import re
import os
import sys
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'checkFASTQ.py',
    description =
    """
    Parses a FASTQ file to check if any sequences don't conform to the SAM spec,
    whether BX:Z: is the last tag in the record, and the counts of: total reads,
    reads without BX:Z: tag, reads with incorrect BX:Z: tag.
    """,
    usage = "checkBAM.py input.bam > output.txt",
    exit_on_error = False
    )

parser.add_argument('input', help = "Input fastq file. Can be gzipped.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

fq_in = args.input

#bxz = re.compile('BX:Z:')
samspec = re.compile('[A-Z][A-Z]:[AifZHB]:')
haplotag = re.compile('A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]')
bxlast = re.compile('BX:Z:A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]$')

with pysam.FastxFile(fq_in) as fh:
    n_reads    = 0
    noBX       = 0
    badBX      = 0
    badSamSpec = 0
    bxNotLast  = 0
    for entry in fh:
        n_reads += 1
        # look for BX:Z: tag
        if 'BX:Z:' in entry.comment:
            # if AXXCXXBXXDXX format isnt found, it's a bad BX
            if not haplotag.search(entry.comment):
                badBX += 1
            splithead = entry.comment.split()
            for i in splithead:
                # if comments dont start with TAG:TYPE:, invalid SAM spec
                if not samspec.match(i):
                    badSamSpec += 1
            # if the BX:Z: isn't at the end, add to bxNotLast
            if splithead[-1].startswith('BX:Z'):
                pass
            else:
                bxNotLast += 1
        else:
            # missing BX:Z: tag
            noBX += 1

values = [str(i) for i in [os.path.basename(fq_in), n_reads, noBX, badBX, badSamSpec, bxNotLast]]
print("\t".join(values), file = sys.stdout)
