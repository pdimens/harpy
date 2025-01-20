#! /usr/bin/env python

import os
import sys
import argparse
from itertools import product
import pysam

parser = argparse.ArgumentParser(
    prog = 'leviathan_bx_shim.py',
    description =
    """
    Deconvolves a BAM file such that BX tags are rewritten to be unique WITHOUT
    hyphenating them. Requires that alignments have already been deconvoled to have
    MI:i tags.
    """,
    usage = "leviathan_bx_shim.py input.bam > output.bam",
    exit_on_error = False
    )

parser.add_argument('input', help = "Input bam file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
if not os.path.exists(args.input):
    parser.error(f"{args.input} was not found")
if not os.path.exists(args.input + ".bai"):
    parser.error(f"{args.input}.bai was not found")
    
# set up a generator for the BX tags
bc_range = [f"{i}".zfill(2) for i in range(1,97)]
bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)
MI_BX = {}

with pysam.AlignmentFile(args.input) as bam_in, pysam.AlignmentFile(sys.stdout.buffer, 'wb', header=bam_in.header) as bam_out:
    # iterate through the bam file
    for record in bam_in.fetch():
        try:
            mi = record.get_tag("MI")
            if mi not in MI_BX:
                try:
                    BX_NEW = "".join(next(bc_generator))          
                except StopIteration:
                    sys.stderr.write(f"Error:\nNumber of unique molecules exceeds the number of possible unique haplotag barcodes ({96**4}).")
                    sys.exit(1)
                MI_BX[mi] = BX_NEW
            record.set_tag("BX", MI_BX[mi])
        except KeyError:
            # no MI tag, just keep the record as-is
            pass
        bam_out.write(record)
