#! /usr/bin/env python
"""Convert inline barcodes into haplotag style ones"""
import os
import sys
import gzip
import argparse
from itertools import zip_longest, product

parser = argparse.ArgumentParser(
    prog = 'inline_to_haplotag.py',
    description = 'Moves inline linked read barcodes to read headers (OX:Z) and converts them into haplotag ACBD format (BX:Z).',
    usage = "inline_to_haplotag.py -f <forward.fq.gz> -r <reverse.fq.gz> -b <barcodes.txt> -p <prefix>",
    exit_on_error = False
    )

parser.add_argument("-f", "--forward", required = True, type = str, help = "Forward reads of paired-end FASTQ file pair (gzipped)")
parser.add_argument("-r", "--reverse", required = True, type = str, help = "Reverse reads of paired-end FASTQ file pair (gzipped)")
parser.add_argument("-l", "--length", required = True, type = int, help = "Length of the barcodes (all must be one length)")
parser.add_argument("-p", "--prefix", required = True, type = str, help = "Prefix for outfile files (e.g. <prefix>.R1.fq.gz)")
parser.add_argument("-b", "--barcodes", required = True, type=str, help="File listing the linked-read barcodes to convert to haplotag format, one barcode per line")
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

def iter_fastq_records(file_handle):
    """Iterate over FASTQ records in a file.
    file_handle: Opened gzip file handle       
    Yields: FASTQ record [header, seq, '+', qual]
    Raises ValueError If file is not in FASTQ format
    """
    record = []
    for line in file_handle:
        line = line.decode().rstrip("\n")
        record.append(line)
        if len(record) == 4:
            # format sanity check
            if not (record[0].startswith("@") and record[2] == "+"):
                raise ValueError("Invalid FASTQ format")
            yield record
            record = []
    if record:
        raise ValueError("Incomplete FASTQ record at end of file")

def validate_barcode(barcode):
    """Validate barcode format (A,C,G,T)."""
    if not set(barcode).issubset({'A','C','G','T'}):
        raise ValueError(f"Invalid barcode format: {barcode}. Barcodes must be captial letters and only contain standard nucleotide values ATCG.")

def process_record(fw_entry, rv_entry, barcode_dict, haplotag_bc, bc_len):
    """convert the barcode to haplotag"""
    # [0] = header, [1] = seq, [2] = +, [3] = qual
    if fw_entry:
        bc_inline = fw_entry[1][:bc_len]
        bc_hap = barcode_dict.get(bc_inline, "A00C00B00D00")
        # the default barcode entry is None, meaning it hasnt been assigned a haplotag equivalent yet
        if not bc_hap:
            bc_hap = "".join(next(haplotag_bc))
            barcode_dict[bc_inline] = bc_hap
        fw_entry[0] = fw_entry[0].split()[0] + f"\tOX:Z:{bc_inline}\tBX:Z:{bc_hap}"
        fw_entry[1] = fw_entry[1][bc_len:]
        fw_entry[3] = fw_entry[3][bc_len:]
        _new_fw = "\n".join(fw_entry) + "\n"
        if rv_entry:
            rv_entry[0]  = rv_entry[0].split()[0] + f"\tOX:Z:{bc_inline}\tBX:Z:{bc_hap}"
            _new_rv = "\n".join(rv_entry) + "\n"
        else:
            _new_rv = None
        return _new_fw, _new_rv
    else:
        # no forward read, therefor no barcode to search for
        rv_entry[0] = rv_entry[0].split()[0] + f"\tBX:Z:A00C00B00D00"
        _new_rv = "\n".join(rv_entry) + "\n"
        return None, _new_rv

bc_range = [f"{i}".zfill(2) for i in range(1,97)]
bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)

bc_dict = {}
# read in barcodes
opener = gzip.open if args.barcodes.lower().endswith('.gz') else open
mode = 'rt' if args.barcodes.lower().endswith('.gz') else 'r'
with opener(args.barcodes, mode) as bc_file:
    for line in bc_file:
        barcode = line.rstrip().split()[0]
        validate_barcode(barcode)
        bc_dict[barcode] = None

# simultaneously iterate the forward and reverse fastq files
fw_out = gzip.open(f"{args.prefix}.R1.fq.gz", "wb", 6)
rv_out = gzip.open(f"{args.prefix}.R2.fq.gz", "wb", 6)

with gzip.open(args.forward, "r") as fw_i, gzip.open(args.reverse, "r") as rv_i:
    for fw_record, rv_record in zip_longest(iter_fastq_records(fw_i), iter_fastq_records(rv_i)):
        new_fw, new_rv = process_record(fw_record, rv_record, bc_dict, bc_generator, args.length)
        if new_fw:
            fw_out.write(new_fw.encode("utf-8"))
        if new_rv:
            rv_out.write(new_rv.encode("utf-8"))

fw_out.close()
rv_out.close()

with open(f"{args.prefix}.barcodes", "w") as bx_out:
    for i,j in bc_dict.items():
        if j:
            bx_out.write(f"{i}\t{j}\n")
