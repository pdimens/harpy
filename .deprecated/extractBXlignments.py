import os
import sys
import pysam

import argparse
parser = argparse.ArgumentParser(
    prog = 'assignMI.py',
    description = 
    """
    Writes all the alignments and their associated molecule to stdout.
    If an alignment does not have a molecule assignment via an MI:i tag,
    it will not be output. This will also ignore alignemnts marked as duplicates.    
    """,
    usage = "extractBXlignments.py INPUT.BAM > OUTPUT",
    exit_on_error = False
    )
parser.add_argument('input', help = "A bam/sam file with MI:i: tags present. If bam, a matching index file should be in the same directory.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

bam_input = args.input

if not os.path.exists(bam_input):
    print(f"Error: {bam_input} not found", file = sys.stderr)
    exit(1)

if bam_input.lower().endswith(".bam"):
    if not os.path.exists(bam_input + ".bai"):
        print(f"Error: {bam_input} requires a matching {bam_input}.bai index file, but one wasn\'t found.", file = sys.stderr)
        exit(1)

with pysam.AlignmentFile(bam_input) as alnfile:
    for record in alnfile.fetch():
        if record.is_unmapped or not record.has_tag("MI") or record.is_duplicate:
            continue
        aln = record.get_blocks()
        if not aln:
            continue

        CHR = record.reference_name
        MI = record.get_tag("MI")
        BX = record.get_tag("BX")
        print(MI, BX, CHR, record.reference_start, record.reference_end, file = sys.stdout)