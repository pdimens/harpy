#! /usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(prog = 'depthWindows.py', description = 'Reads the output of samtools depth -a from stdin and calculates a windowed mean.')
parser.add_argument('windowsize', type= int, help = "The window size to use to calcualte mean depth over (non-overlapping)")

args = parser.parse_args()
_sum = 0
start = 1
end = args.windowsize
lastcontig = None
for line in sys.stdin:
    # Remove the newline character at the end of the line
    line = line.rstrip().split()
    contig = line[0]
    # the contig has changed, make the end position the last position, print output
    if lastcontig and contig != lastcontig:
        winsize = (position + 1) - start
        print(f"{lastcontig}\t{position}\t{_sum / winsize}", file = sys.stdout)
        # reset the window start/end and sum
        _sum = 0
        start = 1
        end = args.windowsize

    position = int(line[1])
    _sum += int(line[2])

    if position == end:
        print(f"{contig}\t{end}\t{_sum / args.windowsize}", file = sys.stdout)
        # reset the window start/end and sum
        _sum = 0
        start = end + 1
        end += args.windowsize
    lastcontig = contig