#! /usr/bin/env python

import os
import re
import sys
import pysam
import random
import argparse

parser = argparse.ArgumentParser(
    prog = 'downsample_bx.py',
    description = 'Downsample a BAM file by values in a record TAG (e.g. BX). Use --invalid-strategy to specify how to handle invalid barcodes: \'keep\' will keep all invalid/missing barcodes without downsampling, \'drop\' will not output any invalid/missing barcodes, \'downsample\' will apply downsampling to invalid/missing barcodes.',
    usage = "downsample_bx.py [-t -s -i] -r 0.3 -o output.bam input.bam",
    )
parser.add_argument('-t', dest = "tag", type=str, default = "BX", help = "SAM tag to use for downsampling (default: %(default)s)")
parser.add_argument('-r', dest = "ratio", type=float, required = True, help = "Downsampling ratio")
parser.add_argument('-i', dest = "strategy", type = str, default = "keep", choices = ["keep","drop", "downsample"], help = "Strategy to handle invalid/missing barcodes (default: %(default)s)")
parser.add_argument('-s', dest = "seed", type = int, help = "Set a random seed (optional)")
parser.add_argument('-o', dest = "outfile", help = "Output bam file")
parser.add_argument('input', type= str, help = "Input bam file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
if not os.path.exists(args.input):
    parser.error(f"Error: {args.input} was not found")
if len(args.tag) != 2:
    parser.error("Error: SAM tags can only be 2 characters long")
if args.ratio < 0 or args.ratio > 1 or args.ratio in [0,1]:
    parser.error("Error: Downsampling ratio must be between 0 and 1")

input_bam = args.input
output_bam = args.outfile
tag_key = args.tag.upper()
downsample_fraction = args.ratio
strat = args.strategy
rng = random.Random(args.seed) if args.seed else random.Random()

#haplotag = re.compile(r'A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}')
invalid = re.compile(r'[AaBbCcDd]00')

# DESIGN
# store records with the same barcode in an array
# then subsample according to the downsample fraction and write to output

with pysam.AlignmentFile(input_bam, "rb") as infile, pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
    record_store_F = []
    record_store_R = []
    last_barcode = None
    for record in infile:
        try:
           barcode = record.get_tag(tag_key)
        except KeyError:
            if strat == "keep":
                outfile.write(record)
            elif strat == "downsample":
                if rng.uniform(0, 1) <= downsample_fraction:
                    outfile.write(record)    
            continue
        if invalid.search(barcode):
            if strat == "keep":
                outfile.write(record)
            elif strat == "downsample":
                if rng.uniform(0, 1) <= downsample_fraction:
                    outfile.write(record)    
            continue
        if not last_barcode or barcode == last_barcode:
            #TODO F/R logic
            record_store.append(record)
            continue
        else:
            for i in record_store:
                if rng.uniform(0, 1) <= downsample_fraction:
                    outfile.write(i)
            # reset the record store
            record_store = []

# TODO ASSESS THAT FORWARD/REVERSE ARE TREATED PROPERLY AS A PAIR
# that is, that if the forward is kept, the reverse is too