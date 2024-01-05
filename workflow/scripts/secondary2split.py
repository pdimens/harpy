#! /usr/bin/env python3

import argparse
import sys
import os
try:
    import pysam
except:
    print("The pysam library is required for operation but was not found, please install it with pip or conda.", file = sys.stderr)
    print("# via pip:", file = sys.stderr)
    print("pip install pysam", file = sys.stderr)
    print("\n# via conda:", file = sys.stderr)
    print("conda install -c bioconda pysam", file = sys.stderr)
    exit(1)

parser = argparse.ArgumentParser(
    prog = 'secondary2split.py',
    description = 'Convert Secondary SAM flags to Split flags for alignments with MAPQ>=30. Writes results to a new BAM file and indexes it with Samtools.',
    usage = "secondary2split.py infile.bam > outfile.bam",
    exit_on_error = False
    )
parser.add_argument("bamfile", help = "Input bam/sam file. If bam format, a corresponding index file (.bam.bai) should be in the same directory.")
parser.add_argument("outfile", help = "Output bam file. A corresponding index file will be created for it.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

if not os.path.exists(args.bamfile):
    print(f"File {args.bamfile} was not found, check the spelling and try again.", file = sys.stderr)
    exit(1)

alnfile = pysam.AlignmentFile(args.bamfile)

with pysam.AlignmentFile(args.outfile, "wb", template = alnfile) as outfile:
    for read in alnfile.fetch():
        if not read.has_tag("XA"):
            if read.mapq >= 30 and read.is_secondary:
                read.is_secondary = False
                read.is_supplementary = True
        _ = outfile.write(read)

alnfile.close()

pysam.index(args.outfile)

exit(0)