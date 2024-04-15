#! /usr/bin/env python3

import re
import sys
import pysam
import argparse

parser = argparse.ArgumentParser(
    prog = 'bxStats.py',
    description = 'Calculate BX molecule length and reads per molecule from BAM file.',
    usage = "bxStats.py input.bam -c cutoff > output.bxstats",
    exit_on_error = False
    )
parser.add_argument('input', help = "Input bam/sam file. If bam, a matching index file should be in the same directory.")
#parser.add_argument('-c','--cutoff', type=int, default = 100000, help = "Distance in base pairs at which alignments with the same barcode should be considered different molecules.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

d = dict()
chromlast = False
alnfile = pysam.AlignmentFile(args.input)

print("barcode_distance", file = sys.stdout)

for read in alnfile.fetch():
    chrm = read.reference_name
    bp   = read.query_alignment_length
    # check if the current chromosome is different from the previous one
    # if so, print the dict to file and empty it (a consideration for RAM usage)
    if chromlast != False and chrm != chromlast:
        # reset the dict
        d = dict()

    if read.is_duplicate or read.is_unmapped:
        continue

    try:
        bx = read.get_tag("BX")
        validBX = False if re.search("[ABCD]0{2,4}", bx) else True
        if not validBX:
            # if invalid bx, skip
            chromlast = chrm
            continue
    except:
        # There is no bx tag, skip
        continue

    aln = read.get_blocks()
    if not aln:
        # unaligned, skip
        continue

    # logic to accommodate split reads 
    # start position of first alignment
    pos_start = aln[0][0]
    # end position of last alignment
    pos_end   = aln[-1][1]

    # create bx entry if it's not present
    if bx not in d.keys():
        d[bx] = {
            "start":  pos_start,
            "end": pos_end,
            "lastpos" : pos_end,
            #"distance" : [],
        }
        chromlast = chrm
        continue

    # distance from last alignment = current aln start - previous aln end
    dist = pos_start - d[bx]["lastpos"]
    if dist > 0:
        print(dist, file = sys.stdout)
    #d[bx]["distance"].append(dist)

    # only if low < currentlow or high > currenthigh
    if pos_start < d[bx]["start"]:
        d[bx]["start"] = pos_start
    if pos_end > d[bx]["end"]:
        d[bx]["end"] = pos_end

    if read.is_reverse or (read.is_forward and not read.is_paired):
        # set the last position to be the end of current alignment
        d[bx]["lastpos"] = pos_end

    # update the chromosome
    chromlast = chrm