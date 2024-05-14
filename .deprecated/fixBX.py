import pysam
import re
import argparse

parser = argparse.ArgumentParser(prog = 'fixBX.py',
                    description = 'Removes any trailing invalid characters incorrectly associated with the BX alignment tag. Example: BX:A01C22B11D84 1:011AAT+TGA will have the invalid text after the full BX tag removed.')
parser.add_argument('input', help = "Input bam/sam file. If bam, a corresponding index file should be in the same directory.")
parser.add_argument('output', help = "Output bam file. This file will also be indexed to create a .bai file.")

args = parser.parse_args()
# the AXXCXXBXXDXX regex match
# haplotag = re.compile("([A-Z]\d{2,}){3,}")

haplotag = re.compile('A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}')

alnfile = pysam.AlignmentFile(args.input)

with pysam.AlignmentFile(args.output, "wb", template = alnfile) as outfile:
    for read in alnfile.fetch():
        try:
            bx = read.get_tag("BX")
            read.set_tag("BX", haplotag.match(bx).group(0))
        except:
            # There is no bx tag
            pass
        outfile.write(read)

alnfile.close()

pysam.index(args.output)

exit(0)