#! /usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(
    prog='makewindows.py',
    description='Create a BED file of fixed intervals from a fasta.fai file (generated with samtools faidx). Nearly identical to bedtools makewindows, except the intervals are nonoverlapping.'
    )
parser.add_argument("-i", dest = "infile", required = True, type=str, metavar = "<input.fasta|fai>", help="input fasta or fasta.fai file")
parser.add_argument("-o", dest = "outfile", required = True, type=str, metavar = "<output.bed>", help="output BED file name")
parser.add_argument("-w", dest = "window", type=int, metavar = "<window_size>", default = 10000, help="interval size (default: %(default)s)")
parser.add_argument("-m", dest = "mode", type=int, metavar = "0 or 1 based", default = 1, help="0 or 1 based intervals (default: %(default)s)")

args = parser.parse_args()
testname = args.infile.lower()
outbed = open(args.outfile, "w")

def readinput(infile, filestream):
    if infile.endswith("gz"):
        return filestream.readline().decode()
    else:
        return filestream.readline()

def makewindows(contig, c_len, windowsize, outfile):
    start = args.mode
    end = windowsize
    starts = [args.mode]
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
            if args.mode == 0:
                c_len -= 1
            makewindows(contig, c_len, args.window, outbed)
        outbed.close()

elif testname.endswith("fasta") or testname.endswith("fa") or testname.endswith("gz"):
    if testname.endswith("gz"):
        import gzip
        fopen = gzip.open(args.infile, "r")
    else:
        fopen = open(args.infile, "r")
    line = readinput(testname, fopen)
    while True:
        c_len=0
        if not line:
            break
        if line.startswith(">"):
            # remove newline, > starting symbol, and any comments after name
            contig = line.rstrip("\n").lstrip(">").split()[0]
            # keep reading until hitting another ">"
            header = False
            while not header:
                line = readinput(testname, fopen)
                if not line:
                    break
                line = line.rstrip("\n")
                if line.startswith(">"):
                    header = True
                    break
                else:
                    c_len += len(line)-1
                    header = False
            if args.mode == 0:
                c_len -= 1
            makewindows(contig, c_len, args.window, outbed)
    outbed.close()
    fopen.close()
else:
    print("Input file not recognized as ending in one of .fai, .fasta, .fa, or .gz (case insensitive).", file = sys.stderr)
    exit(1)
