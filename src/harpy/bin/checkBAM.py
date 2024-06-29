#! /usr/bin/env python

import re
import sys
import os
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'checkBAM.py',
    description =
    """
    Parses an aligment (sam/bam) file to check if the sample name
    matched the RG tag, whether BX:Z: is the last tag in the record,
    and the counts of: total alignments, alignments with an MI:i: tag,
    alignments without BX:Z: tag, incorrect BX:Z: tag.
    """,
    usage = "checkBAM.py input.bam > output.txt",
    exit_on_error = False
    )

parser.add_argument('input', help = "Input bam/sam file. If bam, a matching index file should be in the same directory.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

bam_in = args.input

# regex for EXACTLY AXXCXXBXXDXX
haplotag = re.compile('^A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]$')
bam_pattern = re.compile("\.[bB][aA][mM]$", flags = re.IGNORECASE)

corename = re.sub(bam_pattern, "", os.path.basename(bam_in))

alnfile = pysam.AlignmentFile(bam_in)
if alnfile.header.get("RG")[0]['ID'] == corename:
    nameMismatch = 0
else:
    nameMismatch = 1

n_reads   = 0
noBX      = 0
badBX     = 0
bxNotLast = 0
noMI      = 0

for record in alnfile.fetch():
    n_reads += 1
    tags = [i[0] for i in record.get_tags()]
    # is there a bx tag?
    if 'BX' in tags:
        bx = record.get_tag("BX")
    else:
        noBX += 1
        continue
    # do a regex search to find AXXCXXBXXDXX pattern in the BX
    if not re.search(haplotag, bx):
        # malformed BX tag
        badBX += 1
    # do a search to see if BX:Z: tag is last tag in record
    if tags[-1] != 'BX':
        bxNotLast += 1
    if 'MI' not in tags:
        noMI += 1

alnfile.close()


values = [str(i) for i in [os.path.basename(bam_in), nameMismatch, n_reads, noMI, noBX, badBX, bxNotLast]]
print("\t".join(values), file = sys.stdout)
