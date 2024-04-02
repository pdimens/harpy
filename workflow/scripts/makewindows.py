#! /usr/bin/env python3
import argparse
import sys

parser = argparse.ArgumentParser(
    prog='makewindows.py',
    description='Create a BED file of fixed intervals from a fasta.fai file (generated with samtools faidx). Nearly identical to bedtools makewindows, except the intervals are nonoverlapping.'
    )
parser.add_argument("-i", dest = "infile", required = True, type=str, metavar = "<input.fasta.fai>", help="input fasta.fai file")
parser.add_argument("-o", dest = "outfile", required = True, type=str, metavar = "<output.bed>", help="output BED file name")
parser.add_argument("-w", dest = "window", type=int, metavar = "<window_size>", default = "10000", help="interval size (default: %(default)s)")
args = parser.parse_args()
testname = args.infile.lower()
outbed = open(args.outfile, "w")

def makewindows(contig, c_len, windowsize, outfile):
    start = 1
    end = windowsize
    starts = [1]
    ends = [windowsize]
    while end < c_len:
        end = end + windowsize if (end + windowsize) < c_len else c_len
        ends.append(end)
        start += windowsize
        starts.append(start)
    # write to output file
    for (startpos, endpos) in zip (starts,ends):
        outfile.write(f"{contig}\t{startpos}\t{endpos}\n")

if testname.endswith("fai"):
    with open(args.infile, "r") as fai:
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
            makewindows(contig, c_len, args.window, outbed)
        outbed.close()

elif testname.endswith("fasta") or testname.endswith("fa") or testname.endswith("gz"):
    if testname.endswith("gz"):
        import gzip
        fopen = gzip.open(args.infile, "r")
    else:
        fopen = open(args.infile, "r")
    c_len=0
    contig = ""
    while True:
        # Get next line from file
        if testname.endswith("gz"):
            line = fopen.readline().decode()
        else:
            line = fopen.readline()
        # if line is empty
        # end of file is reached
        if not line:
            break
        if line[0] == ">": 
            if c_len == 0:
                contig = line.rstrip("\n").lstrip(">")
                continue
            else:
                makewindows(contig, c_len, args.window, outbed)
                c_len = 0
                continue
        c_len += len(line)-1
    outbed.close()
    fopen.close()
else:
    print("Input file not recognized as ending in one of .fai, .fasta, .fa, or .gz (case insensitive).", file = sys.stderr)
    exit(1)
