#! /usr/bin/env python3

import pysam
import argparse

parser = argparse.ArgumentParser(prog = 'fixBX.py',
                    description = 'Removes any trailing invalid characters incorrectly associated with the BX alignment tag. Example: BX:A01C22B11D84 1:011AAT+TGA will have the invalid text after the full BX tag removed.')
parser.add_argument('-i', '--input', help = "Input bam/sam file. A corresponding index file should be in the same directory.")
args = parser.parse_args()

alnfile = pysam.AlignmentFile(args.input)
outfile = pysam.AlignmentFile(args.input[0:-4] + ".fix.bam", "wb", template = alnfile)

for read in alnfile.fetch():
    try:
        bx = read.get_tag("BX").split(" ")[0]
        read.set_tag("BX", bx)
    except:
        # There is no bx tag
        pass
    outfile.write(read)

alnfile.close()
outfile.close()
exit(0)