#! /usr/bin/env python
"""Using the ranges of BX molecule start/stop positions, calculates "molecular coverage" across the genome."""

import os
import sys
import gzip
import argparse
from collections import Counter

parser = argparse.ArgumentParser(
    prog = 'molecule_coverage.py',
    description =
    """
    Using the statsfile generated by bx_stats.py from Harpy,
    will calculate "molecular coverage" across the genome.
    Molecular coverage is the "effective" alignment coverage
    if you treat a molecule inferred from linked-read data as
    one contiguous alignment, even though the reads that make
    up that molecule don't cover its entire length. Requires a
    FASTA fai index (like the kind created using samtools faidx)
    to know the actual sizes of the contigs.
    """,
    usage = "molecule_coverage.py -f genome.fasta.fai statsfile > output.cov",
    exit_on_error = False
    )

parser.add_argument('-f', '--fai', required = True, type = str, help = "FASTA index (.fai) file of genome used for alignment")
parser.add_argument('statsfile', help = "stats file produced by harpy via bx_stats.py")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
err = []
for i in [args.statsfile, args.fai]:
    if not os.path.exists(i):
        err.append(i)
if err:
    parser.error("Input files were not found:\n" + ", ".join(err))

# main program
contigs = {}
LASTCONTIG = None

# read the fasta index file as a dict of contig lengths
with open(args.fai, "r", encoding= "utf-8") as fai:
    for line in fai:
        splitline = line.split()
        contig = splitline[0]
        length = splitline[1]
        contigs[contig] = int(length)

def write_coverage(counter_obj, contigname):
    for position,freq in counter_obj.items():
        sys.stdout.write(f"{contigname}\t{position}\t{freq}\n")


with gzip.open(args.statsfile, "rt") as statsfile:
    # read in the header
    line = statsfile.readline()
    # for safety, find out which columns are the contig, start, and end positions
    # just in case this order changes at some point for some reason
    header = line.rstrip().split()
    IDX_CONTIG = None
    IDX_START = None
    IDX_END = None
    for idx,val in enumerate(header):
        if val.strip() == "contig":
            IDX_CONTIG = idx
        if val.strip() == "start":
            IDX_START = idx
        if val.strip() == "end":
            IDX_END = idx
    if IDX_CONTIG is None or IDX_START is None or IDX_END is None:
        sys.stderr.write("Error: Required columns 'contig', 'start', or 'end' not found in header\n")
        sys.exit(1)
    while True:
        line = statsfile.readline()
        if not line:
            if LASTCONTIG:
                # write the last contig to file
                write_coverage(coverage, LASTCONTIG)
            break
        if line.startswith("#"):
            continue
        splitline = line.split()
        contig = splitline[IDX_CONTIG]
        start = int(splitline[IDX_START])
        end = int(splitline[IDX_END])
        if contig != LASTCONTIG:
            if LASTCONTIG:
                # write to file when contig changed
                write_coverage(coverage, LASTCONTIG)
            # reset the counter for the new contig
            coverage = Counter({i: 0 for i in range(1, contigs[contig] + 1)})
        # update the counter with the current row
        coverage.update(range(start, end + 1))
        LASTCONTIG = contig
