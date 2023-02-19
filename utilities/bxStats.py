#! /usr/bin/env python3

import re
import pysam
import argparse

parser = argparse.ArgumentParser(prog = 'bxStats.py',
                    description = 'Calculate BX molecule length and reads per molecule from BAM file.')
parser.add_argument('i', help = "Input bam/sam file. A corresponding index file should be in the same directory.")
args = parser.parse_args()

d = dict()
alnfile = pysam.AlignmentFile(args.i)
outfile = (args.i[0:-4] + ".bx.stats")

for read in alnfile.fetch():
    if read.is_duplicate or read.is_unmapped:
        continue
    try:
        bx = read.get_tag("BX")
        validBX = True
        # do a regex search to find X00 pattern in the BX
        if re.search("[A-Z]0{2,4}", bx):
            # if found, invalid
            bx = "invalidBX"
            validBX = False
    except:
        # There is no bx tag
        bx = "noBX"
        validBX = False
    chrm = read.reference_name
    bp   = read.alen
    if validBX:
        lw, hi = read.blocks[0]
    else:
        lw  = 0
        hi  = 0
    # create chromosome key if not present
    # populate it with the bx stats
    if chrm not in d.keys():
        d[chrm] = {
            bx : {
            "low":  lw,
            "high": hi,
            "bp":   bp,
            "n" :   1
            }
        }
    # create bx stats if it's not present
    elif bx not in d[chrm].keys():
        d[chrm] = {
            bx : {
            "low":  lw,
            "high": hi,
            "bp":   bp,
            "n":    1
            }
        }
    # if BX is present for this chrm, update
    # only if low < currentlow or high > currenthigh
    else:
        if lw < d[chrm][bx]["low"]:
            d[chrm][bx]["low"] = lw
        if hi > d[chrm][bx]["high"]:
            d[chrm][bx]["high"] = hi
        d[chrm][bx]["bp"] += bp
        d[chrm][bx]["n"] += 1


with open(outfile, "w") as fout:
    _ = fout.write("contig\tbx\treads\tstart\tend\tlength_alignment\tlength_inferred\n")
    for chrm in d:
        for bx in d[chrm]:
            _d = d[chrm][bx]
            inferred = str(_d["high"] - _d["low"])
            _ = fout.write(chrm + "\t" + bx + "\t" + str(_d["n"]) + "\t" + str(_d["low"]) + "\t" + str(_d["high"]) + "\t" + str(_d["bp"]) + "\t" + inferred + "\n")