#! /usr/bin/env python
"""Create a BED file of fixed intervals from a fasta or fai file"""
import os
import sys
import gzip
import argparse

parser = argparse.ArgumentParser(
    prog='make_windows.py',
    description='Create a BED file of fixed intervals from a fasta or fai file (generated with samtools faidx). Nearly identical to bedtools makewindows, except the intervals are nonoverlapping.',
    usage = "make_windows.py -w <window.size> -m <0,1> input.fasta > output.bed",
    )
parser.add_argument("input", type=str, help="input fasta or fasta.fai file")
parser.add_argument("-w", "--window", type=int, default = 10000, help="interval size (default: %(default)s)")
parser.add_argument("-m", "--mode", type=int, default = 1, help="0 or 1 based intervals (default: %(default)s)")
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
if not os.path.exists(args.input):
    parser.error(f"{args.input} was not found")
    
testname = args.input.lower()

def is_gzip(file_path):
    """helper function to determine if a file is gzipped, exits if file isn't found"""
    try:
        with gzip.open(file_path, 'rt') as f:
            f.read(10)
        return True
    except gzip.BadGzipFile:
        return False
    except FileNotFoundError:
        sys.stderr.write(f"{file_path} was not found on the system\n")
        sys.exit(1)

def makewindows(_contig, _c_len, windowsize):
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
        sys.stdout.write(f"{_contig}\t{startpos}\t{endpos}\n")

if testname.endswith("fai"):
    with open(args.input, "r", encoding="utf-8") as fai:
        while True:
            line = fai.readline()
            if not line:
                break
            lsplit = line.split("\t")
            contig = lsplit[0]
            c_len = int(lsplit[1])
            c_len += args.mode
            makewindows(contig, c_len, args.window)

elif testname.endswith("fasta") or testname.endswith("fa") or testname.endswith("fasta.gz") or testname.endswith("fa.gz"):
    if is_gzip(args.input):
        fopen = gzip.open(args.input, "rt")
    else:
        fopen = open(args.input, "r", encoding="utf-8")
    line = fopen.readline()
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
                line = fopen.readline()
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
            makewindows(contig, C_LEN, args.window)
    fopen.close()
else:
    sys.stderr.write("Input file not recognized as ending in one of .fai, .fasta[.gz], .fa[.gz] (case insensitive).\n")
    sys.exit(1)
