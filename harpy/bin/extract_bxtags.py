#! /usr/bin/env python

import re
import sys
import pysam
import random
import argparse

parser = argparse.ArgumentParser(
    prog = 'extract_bxtags.py',
    description = 'Extracts all the barcodes present in a SAM/BAM file. Can optionally subsample the barcodes. Use --invalid to specify a proportion of invalid haplotagging to output.',
    usage = "extract_bxtags.py [-i] -b BC -d 15000 input.bam > output.txt",
    )
parser.add_argument('-i','--invalid', metavar = "", type= float, default=1, help = "Proportion of invalid barcodes to include in subsampling/output (default: %(default)s)")
parser.add_argument('-d','--downsample', metavar = "", type=float, help = "Number/fraction of barcodes to downsample to")
parser.add_argument('-r','--random-seed', metavar = "", type= int, help = "Random seed for sampling")
parser.add_argument('-b','--bx-tag', metavar = "", type= str, default="BX", help = "2-character tag with the barcodes (alphanumeric, default: %(default)s)")
parser.add_argument('input', type= str, help = "Input BAM file")
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
if args.invalid < 0 or args.invalid > 1:
    parser.error("--invalid must be between 0 and 1")

bx_tag = args.bx_tag.upper()
inv_prop = args.invalid
rng  = random.Random(args.random_seed) if args.random_seed else random.Random()
invalid_pattern = re.compile(r'[AaBbCcDd]00')
barcodes = set()
mode = "r" if args.input.lower().endswith("sam") else "rb"
with pysam.AlignmentFile(args.input, mode, check_sq=False) as infile:
    for record in infile:
        try:
            barcode = record.get_tag(bx_tag)
            if invalid_pattern.search(barcode):
                # invalid barcode retention
                if rng.random() > inv_prop:
                    continue
        except KeyError:
            continue
        barcodes.add(barcode)

n_bc = len(barcodes)
if args.downsample:
    if args.downsample < 1:
        downsample = int(n_bc * args.downsample)
    else:
        downsample = int(args.downsample)
        if n_bc < downsample:
            raise ValueError(f"The input has fewer barcodes ({n_bc}) than the requested downsampling amount ({args.downsample})")

    for i in rng.sample(sorted(barcodes), downsample):
        sys.stdout.write(f"{i}\n")
else:
    for i in barcodes:
        sys.stdout.write(f"{i}\n")