#! /usr/bin/env python3

import pysam
import re
import argparse

parser = argparse.ArgumentParser(prog = 'filterBXBAM',
                    description = 'Remove alignments from a BAM file that have a least one invalid beadtag barcode.')
parser.add_argument('i', help = "Input bam/sam file. A corresponding index file should be in the same directory.")
args = parser.parse_args()

alnfile = pysam.AlignmentFile(args.i)
outfile = pysam.AlignmentFile(args.i[0:-4] + ".bx.valid.bam", "wb", template = alnfile)

for read in alnfile.fetch():
    try:
        bx = read.get_tag("BX")
    except:
        # There is no bx tag
        continue
    # do a regex search to find X00 pattern in the BX
    if re.search("[A-Z]0{2,4}", bx):
        # if found, invalid and skipped
        continue
    else:
        outfile.write(read)

alnfile.close()
outfile.close()