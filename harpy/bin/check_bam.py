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
    usage = "check_bam.py platform input.bam > output.txt",
    exit_on_error = False
    )

parser.add_argument("platform", metavar='', help= "Linked-read platform\n{10x,haplotagging,stlfr,tellseq}")
parser.add_argument('input', help = "Input bam/sam file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
if not os.path.exists(args.input):
    parser.error(f"{args.input} was not found")
if args.platform not in ["10x","tellseq", "stlfr", "haplotagging"]:
    parser.error("Invalid option for --platform\nMust be one of: 10x, haplotagging, stlfr, tellseq")

if args.platform == "haplotagging":
    bc_pattern = re.compile(r'^A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]')
elif args.platform == "stlfr":
    bc_pattern = re.compile(r'^\d+_\d+_\d+')
else:
    bc_pattern = re.compile(r'^[ATCGN]+')

bam_in = args.input
filename = os.path.basename(bam_in)
bam_pattern = re.compile(r"\.[bB][aA][mM]$", flags = re.IGNORECASE)
corename = re.sub(bam_pattern, "", filename)

N_READS   = 0
NO_BX      = 0
BAD_BX     = 0
BX_NOT_LAST = 0
NO_MI      = 0
NAME_MISMATCH = 0

with pysam.AlignmentFile(bam_in, require_index=False) as alnfile:
    if alnfile.header.get("RG")[0]['ID'] != corename:
        NAME_MISMATCH += 1

    for record in alnfile.fetch(until_eof=True):
        N_READS += 1
        try:
            bx = record.get_tag("BX")
        except KeyError:
            NO_BX += 1
            continue
        # do a regex search to find AXXCXXBXXDXX pattern in the BX
        if not re.search(bc_pattern, bx):
            # malformed BX tag
            BAD_BX += 1
        # do a search to see if BX:Z: tag is last tag in record
        if record.get_tags()[-1][0] != 'BX':
            BX_NOT_LAST += 1
        try:
            mi = record.get_tag("MI")
        except KeyError:
            NO_MI += 1

values = [str(i) for i in [filename, N_READS, NAME_MISMATCH, NO_MI, NO_BX, BX_NOT_LAST, BAD_BX]]
sys.stdout.write("\t".join(values) + "\n")
