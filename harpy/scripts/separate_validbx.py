#! /usr/bin/env python

import argparse
import os
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(
        prog='separate_validbx',
        description='Split a BAM file with BX:Z tags into 2 files, one with valid ACBD barcodes (stdout), one with invalid ACBD barcodes.',
        usage = "separate_validbx invalid.bam input.bam > valid.bam",
        )
    parser.add_argument("invalid_bam", type=str, help="name of output bam for invalid barcodes")
    parser.add_argument("input_bam", type=str, help="input bam file")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not os.path.exists(args.input_bam):
        print(f"{args.input_bam} was provided as input but not found.")
        sys.exit(1)
    if os.path.isdir(args.input_bam):
        print(f"{args.input_bam} is a directory, not a file.")
        sys.exit(1)

    sys.exit(
        subprocess.run(
            ("samtools view -e '[BX]!~\"[ABCD]0{2,4}\"' --unoutput " + f"{args.invalid_bam} {args.input_bam}").split()            
        ).returncode
    )
