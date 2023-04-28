#! /usr/bin/env python3

import pysam
import re
import argparse

parser = argparse.ArgumentParser(prog = 'fixBX.py',
                    description = 'Removes any trailing invalid characters incorrectly associated with the BX alignment tag. Example: BX:A01C22B11D84 1:011AAT+TGA will have the invalid text after the full BX tag removed.')
parser.add_argument('-i', '--input', help = "Input bam/sam file. A corresponding index file should be in the same directory.")
args = parser.parse_args()
# the AXXCXXBXXDXX regex match
haplotag = re.compile('A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}')

alnfile = pysam.AlignmentFile(args.input)
outfile = pysam.AlignmentFile(args.input[0:-4] + ".fix.bam", "wb", template = alnfile)

for read in alnfile.fetch():
    try:
        bx = read.get_tag("BX")
        read.set_tag("BX", haplotag.match(bx).group(0))
    except:
        # There is no bx tag
        pass
    outfile.write(read)

alnfile.close()
outfile.close()
exit(0)