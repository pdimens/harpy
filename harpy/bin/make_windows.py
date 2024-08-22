#! /usr/bin/env python
"""Create a BED file of fixed intervals from a fasta or fai file"""
import argparse
import sys

parser = argparse.ArgumentParser(
    prog='make_windows.py',
    description='Create a BED file of fixed intervals from a fasta or fai file (generated with samtools faidx). Nearly identical to bedtools makewindows, except the intervals are nonoverlapping.'
    )
parser.add_argument("-i", dest = "infile", required = True, type=str, metavar = "<input.fasta|fai>", help="input fasta or fasta.fai file")
parser.add_argument("-o", dest = "outfile", required = True, type=str, metavar = "<output.bed>", help="output BED file name")
parser.add_argument("-w", dest = "window", type=int, metavar = "<window_size>", default = 10000, help="interval size (default: %(default)s)")
parser.add_argument("-m", dest = "mode", type=int, metavar = "0 or 1 based", default = 1, help="0 or 1 based intervals (default: %(default)s)")

args = parser.parse_args()
testname = args.infile.lower()
outbed = open(args.outfile, "w", encoding="utf-8")

def readinput(infile, filestream):
    """automatically read as compressed or not"""
    if infile.endswith("gz"):
        return filestream.readline().decode()
    return filestream.readline()

def makewindows(_contig, _c_len, windowsize, outfile):
    """create a file of the specified windows"""
    start = args.mode
    end = min(_c_len, windowsize)
    starts = [start]
    ends = [end]
    while end < _c_len:
        end = min(end + windowsize, _c_len)
        ends.append(end)
        start += windowsize
        starts.append(start)
    # write to output file
    for (startpos, endpos) in zip (starts,ends):
        outfile.write(f"{_contig}\t{startpos}\t{endpos}\n")

if testname.endswith("fai"):
    with open(args.infile, "r", encoding="utf-8") as fai:
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
            c_len += args.mode
            makewindows(contig, c_len, args.window, outbed)
        outbed.close()

elif testname.endswith("fasta") or testname.endswith("fa") or testname.endswith("gz"):
    if testname.endswith("gz"):
        import gzip
        fopen = gzip.open(args.infile, "r")
    else:
        fopen = open(args.infile, "r", encoding="utf-8")
    line = readinput(testname, fopen)
    while True:
        C_LEN=0
        if not line:
            break
        if line.startswith(">"):
            # remove newline, > starting symbol, and any comments after name
            contig = line.rstrip("\n").lstrip(">").split()[0]
            # keep reading until hitting another ">"
            HEADER = False
            while not HEADER:
                line = readinput(testname, fopen)
                if not line:
                    break
                line = line.rstrip("\n")
                if line.startswith(">"):
                    HEADER = True
                    break
                else:
                    C_LEN += len(line)-1
                    HEADER = False
            C_LEN += args.mode
            makewindows(contig, C_LEN, args.window, outbed)
    outbed.close()
    fopen.close()
else:
    print("Input file not recognized as ending in one of .fai, .fasta, .fa, or .gz (case insensitive).", file = sys.stderr)
    sys.exit(1)
