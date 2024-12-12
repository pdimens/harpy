#! /usr/bin/env python

import os
import re
import sys
import gzip
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'bx_to_end.py',
    description = "Parses a FASTQ or BAM file to move the BX:Z tag to the end of the record.",
    usage = "bx_to_end.py file.[fq|bam]",
    exit_on_error = False
    )

parser.add_argument('input', nargs=1, help = "Input fastq or [indexed] bam file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
infile = args.input[0]
# VALIDATIONS
if not os.path.exists(infile):
    parser.error(f"{infile} was not found")
fq_ext = re.compile(r"\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
if infile.lower().endswith(".bam") or infile.lower().endswith(".sam"):
    is_fastq = False
elif fq_ext.search(infile):
    is_fastq = True
else:
    parser.error(f"Filetype not recognized as one of BAM or FASTQ for file {infile}")

def format_bam(record):
    tags = record.get_tags()
    if "BX" not in tags:
        return record
    bx = record.get_tag("BX")
    # delete the BX tag
    tags = [i for i in tags if i[0] != "BX"]
    # add BX tag to end
    tags.append(("BX", bx))
    record.set_tags(tags)
    return record

def format_fastq(record):
    if "BX:Z" in record.comment:
        splitcomment = record.comment.split()
        # find the BX:Z tag, remove it, and add it to the end
        BX_idx = next((index for index, value in enumerate(splitcomment) if value.startswith("BX:Z")), None)
        bx_tag = splitcomment.pop(BX_idx)
        splitcomment += [bx_tag]
        comment = "\t".join(splitcomment)
    else:
        comment = record.comment
    fastq_req = [
        f"{record.name}\t" + comment,
        record.sequence,
        "+",
        record.quality
    ]
    return "\n".join(fastq_req) + "\n"


if is_fastq:
    with (
        pysam.FastxFile(infile, persist=False) as fq_in,
        gzip.GzipFile(fileobj= sys.stdout.buffer, mode= "wb", compresslevel=6) as fq_out
    ):
        for rec in fq_in:
            fq_out.write(format_fastq(rec).encode())
else:
    try:
        bam_in = pysam.AlignmentFile(infile, "rb")
    except:
        bam_in = pysam.AlignmentFile(infile, "r")
    with pysam.AlignmentFile(sys.stdout.buffer, "wb", template = bam_in) as bam_out:
        for aln_rec in bam_in:
            bam_out.write(format_bam(aln_rec))
    bam_in.close()
