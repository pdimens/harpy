#! /usr/bin/env python
"""parse a fastq file to count BX stats"""
import os
import re
import sys
import argparse
import pysam

def main():
    parser = argparse.ArgumentParser(
        prog = 'count_bx',
        description =
        """
        Parses a FASTQ file to count: total sequences, total number of linked-read barcodes,
        number of valid barcodes, number of invalid BX tags, and a count of positional
        barcode invalidations (e.g. A00, _0_, N)
        """,
        usage = "count_bx platform input.fastq > output.txt",
        exit_on_error = False
        )

    parser.add_argument('platform', choices=["haplotagging", "stlfr", "tellseq"], help = "Type of linked read technology.")
    parser.add_argument('input', help = "Input fastq file. Can be gzipped.")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not os.path.exists(args.input):
        parser.error(f"{args.input} was not found")

    if args.platform == "haplotagging":
        reBARCODE = re.compile(r'BX:Z:(A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2})')
        INVALID_DICT = {
            "A" : 0,
            "C" : 0,
            "B" : 0,
            "D" : 0
        }
        def process_invalid(barcode, segment_dict) -> bool:
            reINVALID = re.compile('([AaBbCcDd])00')
            for i in reINVALID.findall(barcode):
                segment_dict[i] += 1
            return '00' in barcode

    elif args.platform == "stlfr":
        reBARCODE = re.compile('#([0-9]+_[0-9]+_[0-9]+)')
        INVALID_DICT = {
            "1" : 0,
            "2" : 0,
            "3" : 0,
        }
        def process_invalid(barcode, segment_dict):
            split_bc = barcode.split("_")
            for i,j in enumerate(split_bc, 1):
                if j == "0":
                    segment_dict[f'{i}'] += 1
            return '0' in split_bc
    else:
        reBARCODE = re.compile(':([ATCGN]+)')
        INVALID_DICT = dict()
        for i in range(1,19):
            INVALID_DICT[f"{i}"] = 0
        def process_invalid(barcode, segment_dict):
            reINVALID = re.compile('N', flags=re.IGNORECASE)
            for i in reINVALID.finditer(barcode):
                segment_dict[f'{i.start()+1}'] += 1
            return 'N' in barcode

    N_READS = 0
    N_BX = 0
    N_VALID = 0

    with pysam.FastxFile(args.input, persist = False) as fh:
        for entry in fh:
            N_READS += 1
            if args.platform == "haplotagging":
                search_string = entry.comment
            else:
                search_string = entry.name
            
            query = reBARCODE.search(search_string)
            if query:
                N_BX += 1
                barcode = query.group(1)
                is_invalid = process_invalid(barcode, INVALID_DICT)
                if not is_invalid:
                    N_VALID += 1

    sys.stdout.write(f"Reads\t{N_READS}\n")
    sys.stdout.write(f"Barcodes\t{N_BX}\n")
    sys.stdout.write(f"Barcodes_Valid\t{N_VALID}\n")
    sys.stdout.write(f"Barcodes_Invalid\t{N_BX - N_VALID}\n")
    for bx,count in INVALID_DICT.items():
        sys.stdout.write(f"{bx}\t{count}\n")
