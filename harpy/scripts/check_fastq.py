#! /usr/bin/env python
"""parse a fastq file to do BX checks"""
import re
import os
import sys
import argparse
import pysam

class Stats():
    N_READS = 0
    NO_BX = 0
    BAD_BX = 0
    BAD_SAM_SPEC = 0
    BX_NOT_LAST = 0

samspec = re.compile(r'[A-Z][A-Z]:[AifZHB]:')
def check_samspec(fq_comment,_stats):
    splithead = fq_comment.split()
    for i in splithead:
        # if comments dont start with TAG:TYPE:, invalid SAM spec
        if not samspec.match(i):
            _stats.BAD_SAM_SPEC += 1
        # if the BX:Z: isn't at the end, add to BX_NOT_LAST
        if not splithead[-1].startswith('BX:Z'):
            _stats.BX_NOT_LAST += 1

def main():
    parser = argparse.ArgumentParser(
        prog = 'check_fastq',
        description =
        """
        Parses a FASTQ file to check if any sequences don't conform to the SAM spec,
        whether BX:Z: is the last tag in the record, and the counts of: total reads,
        reads without BX:Z: tag, reads with incorrect barcode depending on the platform.
        """,
        usage = "check_fastq platform input.bam > output.txt",
        exit_on_error = False
        )

    parser.add_argument("platform", metavar='', help= "Linked-read platform\n{haplotagging,stlfr,tellseq}")
    parser.add_argument('input', help = "Input fastq file. Can be gzipped.")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not os.path.exists(args.input):
        parser.error(f"{args.input} was not found")
    if args.platform not in ["haplotagging","stlfr","tellseq"]:
        parser.error("Invalid option for --platform\nMust be one of: haplotagging, stlfr, tellseq")

    if args.platform == "haplotagging":
        barcode = re.compile(r'A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]')
        def check_read(fq_record,_stats):
            if 'BX:Z:' not in fq_record.comment:
                _stats.NO_BX += 1
                return
            if not barcode.search(fq_record.comment):
                _stats.BAD_BX += 1
            check_samspec(fq_record.comment,_stats)

    if args.platform == "stlfr":
        barcode = re.compile(r'#\d+_\d+_\d+$')
        def check_read(fq_record,_stats):
            if not barcode.search(fq_record.name):
                _stats.BAD_BX += 1
            check_samspec(fq_record.comment,_stats)

    if args.platform == "tellseq":
        barcode = re.compile(r'\:[ATCGN]+$')
        def check_read(fq_record,_stats):
            if not barcode.search(fq_record.name):
                _stats.BAD_BX += 1
            check_samspec(fq_record.comment,_stats)

    with pysam.FastxFile(args.input, persist=False) as fh:
        for entry in fh:
            Stats.N_READS += 1
            check_read(entry, Stats)

    values = [str(i) for i in [os.path.basename(args.input), Stats.N_READS, Stats.NO_BX, Stats.BAD_BX, Stats.BAD_SAM_SPEC, Stats.BX_NOT_LAST]]
    sys.stdout.write("\t".join(values) + "\n")
