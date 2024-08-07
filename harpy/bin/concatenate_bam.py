#! /usr/bin/env python
"""Concatenate records from haplotagged BAM files"""
import os
import sys
from pathlib import Path
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'concatenate_bam.py',
    description =
    """
    Concatenate records from haplotagged SAM/BAM files while making sure MI:i tags remain unique for every sample.
    This is a means of accomplishing the same as \'samtools cat\', except all MI (Molecule Identifier) tags are updated
    so individuals don't have overlapping MI tags (which would mess up all the linked-read data). You can either provide
    all the files you want to concatenate, or a single file featuring filenames with the \'-b\' option.
    """,
    usage = "concatenate_bam.py -o output.bam file_1.bam file_2.bam...file_N.bam",
    exit_on_error = False
    )
parser.add_argument("alignments", nargs='*', help = "SAM or BAM files")
parser.add_argument("-o", "--out", required = True, type = str, help = "Name of BAM output file")
parser.add_argument("-b", "--bamlist", required = False, type = str, help = "List of SAM or BAM files to concatenate")
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

if (args.alignments and args.bamlist):
    print("Please provide a single file to \'--bamlist\' (-b) featuring all the files you want to concatenate (one per line):", file = sys.stderr)
    print("[example]: concatenate_bam.py -o c_acronotus.bam -b alignments.txt\n", file = sys.stderr)
    
    print("Alternatively, provide the files after \'-o output.bam\':", file = sys.stderr)
    print("[example]: concatenate_bam.py -o c_acronotus.bam sample1.bam sample2.bam", file = sys.stderr)

    sys.exit(1)

if args.bamlist:
    with open(args.bamlist, "r") as bl:
        # read in and filter out commented lines
        aln_list = [i.rstrip() for i in bl.readlines() if not i.startswith("#")]
else:
    if isinstance(args.alignments, str):
        aln_list = [args.alignments]
    else:
        aln_list = args.alignments

# instantiate output file
if aln_list[0].lower().endswith(".bam"):
    if not os.path.exists(f"{aln_list[0]}.bai"):
        pysam.index(aln_list[0])
        # for housekeeping
        DELETE_FIRST_INDEX = True
    else:
        DELETE_FIRST_INDEX = False
with pysam.AlignmentFile(aln_list[0]) as xam_in:
    header = xam_in.header.to_dict()
# Remove all @PG lines
if 'PG' in header:
    del header['PG']
# Add a new @PG line
sys.argv[0] = os.path.basename(sys.argv[0])
new_pg_line = {'ID': 'concatenate', 'PN': 'harpy', 'VN': '1.x', 'CL': " ".join(sys.argv)}
if 'PG' not in header:
    header['PG'] = []
header['PG'].append(new_pg_line)

# update RG lines to match output filename name
header['RG'][0]['ID'] = Path(args.out).stem
header['RG'][0]['SM'] = Path(args.out).stem

with pysam.AlignmentFile(args.out, "wb", header = header) as bam_out:
    # current available unused MI tag
    MI_NEW = 1

    # iterate through the bam files
    for xam in aln_list:
        # create MI dict for this sample
        MI_LOCAL = {}
        # create index if it doesn't exist
        if xam.lower().endswith(".bam"):
            if not os.path.exists(f"{xam}.bai"):
                pysam.index(xam)
                DELETE_INDEX = True
            else:
                DELETE_INDEX = False
        with pysam.AlignmentFile(xam) as xamfile:
            for record in xamfile.fetch():
                try:
                    mi = record.get_tag("MI")
                    # if previously converted for this sample, use that
                    if mi in MI_LOCAL:
                        record.set_tag("MI", MI_LOCAL[mi])
                    else:
                        record.set_tag("MI", MI_NEW)
                        # add to sample conversion dict
                        MI_LOCAL[mi] = MI_NEW
                        # increment to next unique MI
                        MI_NEW += 1
                except:
                    pass
                bam_out.write(record)
        if DELETE_INDEX:
            Path.unlink(f"{xam}.bai")
# just for consistent housekeeping
if DELETE_FIRST_INDEX:
    Path.unlink(f"{aln_list[0]}.bai")
