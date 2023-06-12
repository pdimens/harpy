#! /usr/bin/env python3

import re
import sys
import pysam
import argparse

parser = argparse.ArgumentParser(prog = 'bxStats.py',
                    description = 'Calculate BX molecule length and reads per molecule from BAM file.')
parser.add_argument('input', help = "Input bam/sam file. If bam, a matching index file should be in the same directory.")
parser.add_argument('-c','--cutoff', default = 100000,help = "Distance in base pairs at which alignments with the same barcode should be considered different molecules.")
args = parser.parse_args()

d = dict()
chromlast = False
alnfile = pysam.AlignmentFile(args.input)
#outfile = (args.i[0:-4] + ".bx.stats")

# define write function
# it will only be called when the current alignment's chromosome doesn't
# match the chromosome from the previous alignment
def writestats(x,chr):
    #with open(outf, "a") as fout:
    for bx in x:
        x[bx]["inferred"] = x[bx]["end"] - x[bx]["start"] 
        if x[bx]["mindist"] < 0:
            x[bx]["mindist"] = 0
        outtext = f"{chr}\t{bx}\t" + "\t".join([str(x[bx][i]) for i in ["n", "start","end", "inferred", "bp", "mindist"]])
        #_ = fout.write(outtext + "\n")
        print(outtext, file = sys.stdout)

# write the first line of the output file
#with open(outfile, "w") as fout:
#    _ = fout.write("contig\tbx\treads\tstart\tend\tlength_inferred\taligned_bp\tmindist\texceed_50k\n")
print("contig\tbx\treads\tstart\tend\tlength_inferred\taligned_bp\tmindist", file = sys.stdout)

#alnfile = pysam.AlignmentFile("/home/pdimens/subset.bam")


for read in alnfile.fetch():
    chrm = read.reference_name
    bp   = read.query_alignment_length
    # check if the current chromosome is different from the previous one
    # if so, print the dict to file and empty it (a consideration for RAM usage)
    if chromlast != False and chrm != chromlast:
        writestats(d, chromlast)
        d = dict()
    if read.is_duplicate or read.is_unmapped:
        continue
    try:
        bx = read.get_tag("BX")
        validBX = True
        # do a regex search to find X00 pattern in the BX
        if re.search("[ABCD]0{2,4}", bx):
            # if found, invalid
            bx = "invalidBX"
            validBX = False
    except:
        # There is no bx tag
        bx = "noBX"
        validBX = False
    if validBX:
        pos_start, pos_end = read.blocks[0]
    else:
        pos_start  = 0
        pos_end  = 0

    # create bx stats if it's not present
    if bx not in d.keys():
        d[bx] = {
            "start":  pos_start,
            "end": pos_end,
            "bp":   bp,
            "n":    1,
            "lastpos" : pos_end,
            "mindist" : -1,
            "current_suffix": 0
        }
        chromlast = chrm
        continue

    # if invalid/absent BX, skip the distance stuff
    if bx in ["noBX", "invalidBX"]:
        chromlast = chrm
        continue

    orig = bx
    # section to add current suffix to barcode name, if exists
    if d[orig]["current_suffix"] > 0:
        bx = orig + "." + str(d[orig]["current_suffix"])

    # distance from last alignment = current aln start - previous aln end
    if read.is_forward or (read.is_reverse and not read.is_paired):
        dist = pos_start - d[bx]["lastpos"]
        # if the distance between alignments is >100kbp, it's a different molecule
        # so we'll +1 the suffix of the main barcode and relabel this one as BX + suffix
        # since it's a new entry, we initialize it and move on
        if dist > args.cutoff:
            d[orig]["current_suffix"] += 1
            bx = orig + "." + str(d[orig]["current_suffix"])
            d[bx] = {
                "start":  pos_start,
                "end": pos_end,
                "bp":   bp,
                "n":    1,
                "lastpos" : pos_end,
                "mindist" : -1,
                "current_suffix": 0
            }
            chromlast = chrm
            continue 
        if dist < d[bx]["mindist"] or d[bx]["mindist"] < 0:
            d[bx]["mindist"] = dist

    # if BX is present, update
    d[bx]["bp"] += bp
    d[bx]["n"]  += 1
    chromlast = chrm

    # only if low < currentlow or high > currenthigh
    if pos_start < d[bx]["start"]:
        d[bx]["start"] = pos_start
    if pos_end > d[bx]["end"]:
        d[bx]["end"] = pos_end

    if read.is_reverse or (read.is_forward and not read.is_paired):
        # set the last position to be the end of current alignment
        d[bx]["lastpos"] = pos_end



