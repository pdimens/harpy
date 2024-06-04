#! /usr/bin/env python3
import argparse
parser = argparse.ArgumentParser(
    prog='makewindows.py',
    description='Create a BED file of fixed intervals from a fasta.fai file (generated with samtools faidx). Nearly identical to bedtools makewindows, except the intervals are nonoverlapping.'
    )
parser.add_argument("-i", dest = "infile", required = True, type=str, metavar = "<input.fasta.fai>", help="input fasta.fai file")
parser.add_argument("-o", dest = "outfile", required = True, type=str, metavar = "<output.bed>", help="output BED file name")
parser.add_argument("-w", dest = "window", type=int, metavar = "<window_size>", default = "10000", help="interval size (default: %(default)s)")
args = parser.parse_args()

with open(args.infile, "r") as fai:
    outbed = open(args.outfile, "w")
    while True:
        # Get next line from file
        line = fai.readline()
        # if line is empty
        # end of file is reached
        if not line:
            break
        # split the line by tabs
        lsplit = line.split("\t")
        contig = lsplit[0]
        c_len = int(lsplit[1])
        start = 1
        end = args.window
        starts = [1]
        ends = [args.window]
        while end < c_len:
            end = end + args.window if (end + args.window) < c_len else c_len
            ends.append(end)
            start += args.window
            starts.append(start)

        # write to output file
        for (startpos, endpos) in zip (starts,ends):
            outbed.write(f"{contig}\t{startpos}\t{endpos}\n")
    outbed.close()
