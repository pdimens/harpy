#! /usr/bin/env python
"""Parse a phase block file from HapCut2 to pull out summary information"""

import os
import sys
import argparse

parser = argparse.ArgumentParser(
    prog='parse_phaseblocks.py',
    description='Parse a phase block file from HapCut2 to pull out summary information',
    usage = "parse_phaseblocks.py input > output.txt"
    )
parser.add_argument("input", type=str, help="input HapCut2 phase blocks file")
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
if not os.path.exists(args.input):
    parser.error(f"{args.input} was not found\n")

samplename = args.input.replace(".blocks", "")

with open(args.input, "r", encoding="utf-8") as blocks:
    FIRST_LINE = True
    for line in blocks:
        lsplit = line.split()
        if lsplit[0] == "BLOCK:":
            n_snp     = int(lsplit[6])
            len_block = int(lsplit[8])
            FIRST_LINE = True
        else:
            if FIRST_LINE:
                pos_start = int(lsplit[4])
                contig    = lsplit[3]
                FIRST_LINE = False
                sys.stdout.write(f"{samplename}\t{contig}\t{n_snp}\t{pos_start}\t{len_block}\n")
            else:
                continue
