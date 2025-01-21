#! /usr/bin/env python
"""Concatenate records from haplotagged BAM files"""
import os
import sys
from pathlib import Path
import argparse
from itertools import product
import pysam

parser = argparse.ArgumentParser(
    prog = 'concatenate_bam.py',
    description =
    """
    Concatenate records from haplotagged SAM/BAM files while making sure MI:i tags remain unique for every sample.
    This is a means of accomplishing the same as \'samtools cat\', except all MI (Molecule Identifier) tags are updated
    so individuals don't have overlapping MI tags (which would mess up all the linked-read info). You can either provide
    all the files you want to concatenate, or a single file featuring filenames with the \'-b\' option. Use the \'--bx\'
    option to also rewrite BX tags such that they are unique between samples too.
    """,
    usage = "concatenate_bam.py [--bx] file_1.bam..file_N.bam > out.bam",
    exit_on_error = False
    )
parser.add_argument("alignments", nargs='*', help = "SAM or BAM files")
#parser.add_argument("-o", "--out", required = True, type = str, help = "Name of BAM output file")
parser.add_argument("-b", "--bamlist", required = False, type = str, help = "List of SAM or BAM files to concatenate")
parser.add_argument("--bx", dest = "bx_unique", action='store_true', help="Also rewrite BX tags for uniqueness")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

if (args.alignments and args.bamlist):
    sys.stderr.write("Please provide a single file to \'--bamlist\' (-b) featuring all the files you want to concatenate (one per line):\n")
    sys.stderr.write("[example]: concatenate_bam.py -b alignments.txt > c_acronotus.bam\n\n")
    sys.stderr.write("Alternatively, provide the files at the end of the command:\n",)
    sys.stderr.write("[example]: concatenate_bam.py sample1.bam sample2.bam > c_acronotus.bam\n")
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

# validate files exist
err = []
for i in aln_list:
    if not os.path.exists(i):
        err.append(i)
if err:
    parser.error("Some input files were not found on the system:\n" + ", ".join(err))

# Get the max number of unique haplotag barcodes
haplotag_limit = 96**4

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
header['RG'][0]['ID'] = "concat"
header['RG'][0]['SM'] = "concat"

# set up a generator for the BX tags if --bx was invoked
if args.bx_unique:
    bc_range = [f"{i}".zfill(2) for i in range(1,97)]
    bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)

with pysam.AlignmentFile(sys.stdout.buffer, "wb", header = header) as bam_out:
    # current available unused MI tag
    MI_NEW = 1
    if args.bx_unique:
        BX_NEW = "".join(next(bc_generator))
    else:
        BX_NEW = None
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
                        record.set_tag("MI", MI_LOCAL[mi][0])
                        if args.bx_unique:
                            record.set_tag("BX", MI_LOCAL[mi][1])
                    else:
                        record.set_tag("MI", MI_NEW)
                        if args.bx_unique:
                            record.set_tag("BX", BX_NEW)
                        # add to sample conversion dict
                        MI_LOCAL[mi] = [MI_NEW, BX_NEW]
                        # increment to next unique MI
                        MI_NEW += 1
                        if args.bx_unique:
                            if MI_NEW > haplotag_limit:
                                raise IndexError
                            BX_NEW = "".join(next(bc_generator))
                except IndexError:
                    errtext = f"Error:\nNumber of unique molecules exceeds the number of possible unique haplotag barcodes (96^4 = {haplotag_limit}). "
                    errtext += "Consider pooling fewer individuals per group.\n\nIf this concatenation was performed as part of the harpy sv leviathan workflow, "
                    errtext += "there is a limitation in LEVIATHAN where it does not recognize hyphenated (deconvolved) linked-read barcodes, which necessitates using all possible unique standard haplotag barcodes. "
                    errtext += "Consider using the 'sv naibr' workflow, which uses unique MI tags instead.\n"
                    sys.stderr.write(errtext)
                    sys.exit(1)
                except KeyError:
                    pass
                bam_out.write(record)
        if DELETE_INDEX:
            Path.unlink(f"{xam}.bai")
# just for consistent housekeeping
if DELETE_FIRST_INDEX:
    Path.unlink(f"{aln_list[0]}.bai")
