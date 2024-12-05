#! /usr/bin/env python

import os
import re
import sys
import pysam
import random
import argparse

desc = """Downsample a BAM file by values in a record TAG (e.g. BX).
== Downsampling Methods ==
- value of -d <1: downsamples reads within a barcode
    - drops unpaired reads
    - e.g. '-d 0.3' downsamples reads sharing a single barcode by 0.3

- value of -d >1: downsamples the barcodes themselves
    - retains unpaired reads
    - e.g. '-d 1000' retains all reads associated with 1000 random barcodes

== Invalid/Missing Barcodes ==
Use -i to specify how to handle invalid barcodes:
- 'keep': keep all invalid/missing barcodes
- 'drop': don't output any invalid/missing barcodes
- 'downsample': (-d <1 only) also downsample invalid/missing barcodes
"""

parser = argparse.ArgumentParser(
    prog = 'downsample_bx.py',
    description = desc,
    usage = "downsample_bx.py [-t -s -i] -d DOWNSAMPLE -o OUTFILE input.bam",
    formatter_class  =  argparse.RawDescriptionHelpFormatter
    )
parser.add_argument('-d', dest = "downsample", type=float, required = True, help = "Downsampling amount")
parser.add_argument('-i', dest = "invalid", type = str, default = "keep", choices = ["keep","drop", "downsample"], help = "Strategy to handle invalid/missing barcodes (default: %(default)s)")
parser.add_argument('-s', dest = "seed", type = int, help = "Set a random seed (optional)")
parser.add_argument('-t', dest = "tag", type=str, default = "BX", help = "SAM tag to use for downsampling (default: %(default)s)")
parser.add_argument('-o', dest = "outfile", help = "Output bam file")
parser.add_argument('input', type= str, help = "Input bam file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
if not os.path.exists(args.input):
    parser.error(f"{args.input} was not found")
if len(args.tag) != 2:
    parser.error("SAM tags can only be 2 characters long")
if args.downsample <= 0:
    parser.error("Downsampling ratio must be greater than 0")
elif args.downsample >= 1:
    mode = "across"
else:
    mode = "within"

input_bam = args.input
output_bam = args.outfile
tag_key = args.tag.upper()
downsample_fraction = args.downsample if args.downsample < 1 else int(args.downsample)
strat = args.invalid
rng = random.Random(args.seed) if args.seed else random.Random()

#haplotag = re.compile(r'A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}')
invalid = re.compile(r'[AaBbCcDd]00')

# DESIGN
# store records with the same barcode in an array
# then subsample according to the downsample fraction and write to output
if mode == "within":
    with pysam.AlignmentFile(input_bam, "rb") as infile, pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
        record_store_F = []
        record_store_R = []
        last_barcode = None
        for record in infile:
            try:
                barcode = record.get_tag(tag_key)
                if isinstance(barcode, int):
                    pass
                elif invalid.search(barcode):
                    if strat == "keep":
                        outfile.write(record)
                    elif strat == "downsample":
                        if rng.uniform(0, 1) <= downsample_fraction:
                            outfile.write(record)    
                    continue
            except KeyError:
                if strat == "keep":
                    outfile.write(record)
                elif strat == "downsample":
                    if rng.uniform(0, 1) <= downsample_fraction:
                        outfile.write(record)    
                continue
            if not last_barcode or barcode == last_barcode:
                if record.is_forward and record.is_paired:
                    record_store_F.append(record)
                elif record.is_reverse and record.is_paired:
                    record_store_R.append(record)         
                continue
            else:
                for i,j in zip(record_store_F, record_store_R):
                    if rng.uniform(0, 1) <= downsample_fraction:
                        outfile.write(i)
                        outfile.write(j)
                # reset the record stores
                record_store_F = []
                record_store_R = []
else:
    # Design: read input file, get list of valid barcodes, subsample valid barcodes, read input again, output only sampled barcodes
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        barcodes = set()
        for record in infile:
            try:
                barcode = record.get_tag(tag_key)
                if isinstance(barcode, int):
                    pass # an int from an MI-type tag
                elif invalid.search(barcode):
                    continue
            except KeyError:
                continue
            barcodes.add(barcode)
        n_bc = len(barcodes)
        if downsample_fraction > n_bc:
            sys.stderr.write(f"The number of intended barcodes to downsample to ({downsample_fraction}) is greater than the number of barcodes in the input file ({n_bc}).\n")
            sys.exit(1)
        barcodes = random.sample(sorted(barcodes), downsample_fraction)
    with pysam.AlignmentFile(input_bam, "rb") as infile, pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
        for record in infile:
            try:
                barcode = record.get_tag(tag_key)
                if isinstance(barcode, int):
                    pass # an int from an MI-type tag
                elif invalid.search(barcode):
                    if strat == "keep":
                        outfile.write(record)
                    continue
            except KeyError:
                if strat == "keep":
                    outfile.write(record)
                continue
            if barcode in barcodes:
                outfile.write(record)
