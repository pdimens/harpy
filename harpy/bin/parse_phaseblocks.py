#! /usr/bin/env python
"""Parse a phase block file from HapCut2 to pull out summary information"""
import argparse
import sys

parser = argparse.ArgumentParser(
    prog='parse_phaseblocks.py',
    description='Parse a phase block file from HapCut2 to pull out summary information'
    )
parser.add_argument("-i", dest = "infile", required = True, type=str, metavar = "<file.blocks>", help="input HapCut2 phase blocks file")
args = parser.parse_args()

samplename = args.infile.replace(".blocks", "")

with open(args.infile, "r", encoding="utf-8") as blocks:
    FIRST_LINE = True
    while True:
        # Get next line from file
        line = blocks.readline()
        # if line is empty, end of file is reached
        if not line:
            break

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
                print(f"{samplename}\t{contig}\t{n_snp}\t{pos_start}\t{len_block}", file = sys.stdout)
            else:
                continue
