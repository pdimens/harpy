#! /usr/bin/env python
"""Convert 10X style barcodes into Haplotag style ones"""
import os
import sys
import gzip
import argparse
from itertools import zip_longest, product

parser = argparse.ArgumentParser(
    prog = '10xtoHaplotag.py',
    description = 'Converts 10x linked reads to haplotag linked reads with barcodes in BX:Z: and OX:Z: header tags.',
    usage = "10xtoHaplotag.py -f <forward.fq.gz> -r <reverse.fq.gz> -b <barcodes.txt> -p <prefix> > barcodes.conversion.txt",
    exit_on_error = False
    )

parser.add_argument("-f", "--forward", required = True, type = str, help = "Forward reads of paired-end FASTQ file pair (gzipped)")
parser.add_argument("-r", "--reverse", required = True, type = str, help = "Reverse reads of paired-end FASTQ file pair (gzipped)")
parser.add_argument("-p", "--prefix", required = True, type = str, help = "Prefix for outfile FASTQ files (e.g. <prefix>.R1.fq.gz)")
parser.add_argument("-b", "--barcodes", required = True, type=str, help="File listing the 10X barcodes to convert to haplotag format, one barcode per line")
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
err = []
for i in [args.forward, args.reverse, args.barcodes]:
    if not os.path.exists(i):
        err.append(i)
if err:
    parser.error("Some input files were not found on the system:\n" + ", ".join(err))

def process_record(fw_entry, rv_entry):
    """convert the 10X to haplotag"""
    # [0] = header, [1] = seq, [2] = +, [3] = qual
    bc10x = fw_entry[1][:16]
    bchap = bc_dict.get(bc10x, "A00C00B00D00")
    if not bchap:
        bchap = "".join(next(bc_generator))
        bc_dict[bc10x] = bchap
    _new_fw  = fw_entry[0].split()[0] + f"\tOX:Z:{bc10x}\tBX:Z:{bchap}\n"
    _new_fw += fw_entry[1][16:] + "\n"
    _new_fw += fw_entry[2] + "\n"
    _new_fw += fw_entry[3][16:] + "\n"
    _new_rv  = rv_entry[0].split()[0] + f"\tOX:Z:{bc10x}\tBX:Z:{bchap}\n"
    _new_rv += "\n".join(rv_entry[1:3])
    return _new_fw, _new_rv

bc_range = [f"{i}".zfill(2) for i in range(1,97)]
bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)

bc_dict = {}

# read in barcodes
with open(args.barcodes, "r") as bc_file:
    while True:
        # Get next line from file
        if args.barcodes.endswith("gz"):
            line = bc_file.readline().decode()
        else:
            line = bc_file.readline()
        # if line is empty
        # end of file is reached
        if not line:
            break
        bc = line.rstrip("\n").split()
        _10X = str(bc[0])
        bc_dict[_10X] = None

# simultaneously iterate the forward and reverse fastq files
fw_reads = args.forward
rv_reads = args.reverse

fw_out = gzip.open(f"{args.prefix}.R1.fq.gz", "wb", 6)
rv_out = gzip.open(f"{args.prefix}.R2.fq.gz", "wb", 6)

with gzip.open(fw_reads, "r") as fw_i, gzip.open(rv_reads, "r") as rv_i:
    # store the fq records here
    fw_record = []
    rv_record = []
    i = 0
    for fw, rv in zip_longest(fw_i, rv_i):
        fw_line = fw.decode().rstrip("\n")
        rv_line = rv.decode().rstrip("\n")
        if fw_line.startswith("@") and i > 0:
            i += 1
            # process the full record
            new_fw, new_rv = process_record(fw_record, rv_record)
            # write new record to files
            fw_out.write(new_fw.encode("utf-8"))
            rv_out.write(new_rv.encode("utf-8"))
            # reset the record
            fw_record = [fw_line]
            rv_record = [rv_line]
        elif fw_line.startswith("@") and i == 0:
            # its the first record don't write anything
            i += 1
            fw_record = [fw_line]
            rv_record = [rv_line]
        else:
            # just append the other lines to the record
            fw_record.append(fw_line)
            rv_record.append(rv_line)

fw_out.close()
rv_out.close()

for i,j in bc_dict.items():
    if j:
        print(f"{i}\t{j}", file = sys.stdout)
