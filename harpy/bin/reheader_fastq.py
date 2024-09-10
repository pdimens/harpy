#! /usr/bin/env python

import re
import os
import sys
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'reheader_fastq.py',
    description =
    """
    Reformat the sequence headers of NGSNGS output to match the format of
    DWGSIM. For use in the harpy simulate linkedreads workflow.
    """,
    usage = "reheader_fastq.py input.fq > output.fq",
    exit_on_error = False
    )

parser.add_argument('input', help = "Input fastq file. Can be gzipped.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

with pysam.FastxFile(args.input, persist=False) as fh:
    CONVERSION = [2,3,3,4,5] # etc...
    for entry in fh:
        ngsngs_split = re.split(r"[_\:]", entry.name)
        dwgsim_split = [ngsngs_split[i] for i in CONVERSION]
        dwgsim_name = # MANUAL JOIN HERE
        print(dwgsim_name, file = sys.stdout)
        print(entry.sequence,file = sys.stdout)
        print(entry.quality, file = sys.stdout)
        print("+", file = sys.stdout)