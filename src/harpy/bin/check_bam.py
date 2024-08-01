#! /usr/bin/env python
"""parse a bam file to check for BX stats"""
import re
import sys
import os
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'check_bam.py',
    description =
    """
    Parses an aligment (sam/bam) file to check if the sample name
    matched the RG tag, whether BX:Z: is the last tag in the record,
    and the counts of: total alignments, alignments with an MI:i: tag,
    alignments without BX:Z: tag, incorrect BX:Z: tag.
    """,
    usage = "check_bam.py input.bam > output.txt",
    exit_on_error = False
    )

parser.add_argument('input', help = "Input bam/sam file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

bam_in = args.input
if bam_in.lower().endswith(".bam"):
    if not os.path.exists(bam_in + ".bai"):
        pysam.index(bam_in)

# regex for EXACTLY AXXCXXBXXDXX
haplotag = re.compile('^A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]')
bam_pattern = re.compile(r"\.[bB][aA][mM]$", flags = re.IGNORECASE)

corename = re.sub(bam_pattern, "", os.path.basename(bam_in))

alnfile = pysam.AlignmentFile(bam_in)
if alnfile.header.get("RG")[0]['ID'] == corename:
    NAME_MISMATCH = 0
else:
    NAME_MISMATCH = 1

N_READS   = 0
NO_BX      = 0
BAD_BX     = 0
BX_NOT_LAST = 0
NO_MI      = 0

for record in alnfile.fetch():
    N_READS += 1
    tags = [i[0] for i in record.get_tags()]
    # is there a bx tag?
    if 'BX' in tags:
        bx = record.get_tag("BX")
    else:
        NO_BX += 1
        continue
    # do a regex search to find AXXCXXBXXDXX pattern in the BX
    if not re.search(haplotag, bx):
        # malformed BX tag
        BAD_BX += 1
    # do a search to see if BX:Z: tag is last tag in record
    if tags[-1] != 'BX':
        BX_NOT_LAST += 1
    if 'MI' not in tags:
        NO_MI += 1

alnfile.close()


values = [str(i) for i in [os.path.basename(bam_in), NAME_MISMATCH, N_READS, NO_MI, NO_BX, BAD_BX, BX_NOT_LAST]]
print("\t".join(values), file = sys.stdout)
