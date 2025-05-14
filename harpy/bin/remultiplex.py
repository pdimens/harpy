#! /usr/bin/env python

import argparse
from itertools import product
import os
import re
import subprocess
import sys
import pysam

parser = argparse.ArgumentParser(
    prog = 'remultiplex.py',
    description =
    """
    Converts the linked-read barcode into nucleotide format (if necessary) and adds it to the beginning
    of the sequence. Writes a file of the barcode conversion map.
    """,
    usage = "remultiplex.py prefix input.r1.fq input.r2.fq",
    exit_on_error = False
    )

parser.add_argument("-s", "--scan", default=100, type = int, help = "Number of reads to scan to identify barcode location and format")
parser.add_argument("-i", "--preserve-invalid", type = int, help = "Retain the uniqueness of invalid barcodes")
parser.add_argument("-m", "--map", type = int, help = "Write a map for the barcode-to-nucleotide conversion")
parser.add_argument("prefix", help= "Filename prefix for output")
parser.add_argument('r1_fq',  help = "Input fastq R1 file. Can be gzipped.")
parser.add_argument('r2_fq', help = "Input fastq R2 file. Can be gzipped.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

def bx_barcode(rec):
    bx = [i for i in rec.comment.split() if i.startswith("BX:Z:")]
    #bx = re.search(r"BX:Z:[^\s]*(?=\s)", rec.comment)
    if bx:
        return bx[0].removeprefix("BX:Z:")
    else:
        return None

def inline_barcode(rec):
    bx = re.search(r"(?:\:[ATCGN]+$|#\d+_\d+_\d+$)", rec.name)
    if bx:
        return bx[0][1:]
    else:
        return None

def is_invalid(bx):
    return "0" in bx.split("_") or bool(re.search(r"(?:N|[ABCD]00)", bx))

NUCLEOTIDE_FMT = False
## find the barcode format
with pysam.FastqFile(args.r1_fq, "r") as in_fq:
    for n,record in enumerate(in_fq, 1):
        if bx_barcode(record):
            bx_search = bx_barcode
            _bx = bx_search(record)
            if re.search(r"^[ATCGN]+$", _bx):
                NUCLEOTIDE_FMT = True
            break
        elif inline_barcode(record):
            bx_search = inline_barcode
            _bx = bx_search(record)
            if re.search(r"^[ATCGN]+$", _bx):
                NUCLEOTIDE_FMT = True
            break
        if n > args.scan:
            print(f"Scanned the first {args.scan} reads of {os.path.basename(args.r1_fq)} and was unable to locate barcodes in the BX:Z field nor as a TELLseq or stLFR suffix in the read ID.")
            sys.exit(1)

bc_inventory = {}
bc_iter = product(*["ATCG" for i in range(18)])
bc_iter_inv = product(*(["N"] + ["ATCGN" for i in range(17)]))

for i,fq in enumerate([args.r1_fq, args.r2_fq],1):
    with pysam.FastqFile(fq, "r") as in_fq, open(f"{args.prefix}.R{i}.fq.gz", "wb") as out_fq:
        gzip = subprocess.Popen(["gzip"], stdin = subprocess.PIPE, stdout = out_fq)
        for record in in_fq:
            _bx = bx_search(record)
            if not _bx:
                record.sequence = "N"*18 + "NNNN" + record.sequence
                record.quality  = "!"*16 + "!!!!" + record.quality
            else:
                if NUCLEOTIDE_FMT:
                    record.sequence = _bx + "NNNN" + record.sequence
                    record.quality =  "I"*len(_bx) + "!!!!" + record.quality
                else:
                    nuc_bx = bc_inventory.get(_bx, None)
                    if not nuc_bx:
                        if is_invalid(_bx):
                            nuc_bx = "".join(next(bc_iter_inv)) if args.preserve_invalid else "N"*18
                        else:
                            nuc_bx = "".join(next(bc_iter))
                        bc_inventory[_bx] = nuc_bx
                    record.sequence = nuc_bx + "NNNN" + record.sequence
                    record.quality =  "I"*18 + "!!!!" + record.quality
            gzip.stdin.write(str(record).encode("utf-8") + b"\n")

if args.map:
    with open(f"{prefix}.barcode.map", "w") as bc_out:
        for bx,nuc in bc_inventory.items():
            bc_out.write(f"{bx}\t{nuc}\n")