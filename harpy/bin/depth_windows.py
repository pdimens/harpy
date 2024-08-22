#! /usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(
    prog = 'depth_windows.py',
    description = 'Reads the output of samtools depth -a from stdin and calculates a windowed mean',
    usage = "samtools depth -a file.bam | depth_windows.py windowsize > output.txt",
    )
parser.add_argument('windowsize', type= int, help = "The window size to calcualte mean depth over (non-overlapping)")
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
if args.windowsize < 1:
    parser.error("Error: window size must be greater than 0")
if args.windowsize == 1:
    # just print the input to output
    for line in sys.stdin:
        sys.stdout.write(line)
    sys.exit(0)


_SUM = 0
START = 1
END = args.windowsize
LAST_CONTIG = None
POSITION = 0

for line in sys.stdin:
    # Remove the newline character at the END of the line
    line = line.rstrip().split()
    contig = line[0]
    # the contig has changed, make the END POSITION the last POSITION, print output
    if LAST_CONTIG and contig != LAST_CONTIG:
        WINSIZE = (POSITION + 1) - START
        if WINSIZE > 0:
            depth = _SUM / WINSIZE
            sys.stdout.write(f"{LAST_CONTIG}\t{POSITION}\t{depth}\n")
        # reset the window START/END and sum
        _SUM = 0
        START = 1
        END = args.windowsize

    POSITION = int(line[1])
    _SUM += int(line[2])

    if POSITION == END:
        depth = _SUM / args.windowsize
        sys.stdout.write(f"{contig}\t{END}\t{depth}\n")
        # reset the window START/END and sum
        _SUM = 0
        START = END + 1
        END += args.windowsize
    LAST_CONTIG = contig
