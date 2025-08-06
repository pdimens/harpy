#! /usr/bin/env python
"""Rename a sam/bam file and modify the @RG tag of the alignment file to reflect the change for both ID and SM."""
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(
        prog='rename_bam',
        description='Rename a sam/bam file and modify the @RG tag of the alignment file to reflect the change for both ID and SM. This process creates a new file \'newname.bam\' and you may use -d to delete the original file. Requires samtools.',
        usage = "rename_bam [-d] new_name input.bam"
        )
    parser.add_argument("name", type=str, help="new sample name")
    parser.add_argument("input", type=str, help="input bam or sam file")
    parser.add_argument("-d", "--delete", dest = "delete", action='store_true', help="delete the original file")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if not os.path.exists(args.input):
        parser.error(f"{args.input} was not found")

    outdir = os.path.dirname(args.input)

    JOB_STATUS = os.system(
        f"samtools addreplacerg -r \"ID:{args.name}\\tSM:{args.name}\" -o {outdir}/{args.name}.bam {args.input}"
        )

    if JOB_STATUS != 0:
        sys.stderr.write("samtools addreplacerg failed with an error\n")  
        sys.exit(JOB_STATUS)
    else:
        if args.delete:
            try:
                os.remove(args.input)
            except OSError:
                sys.stderr.write(f"Failed to delete {args.input}, but otherwise samtools was successful.\n")
        sys.exit(JOB_STATUS)
