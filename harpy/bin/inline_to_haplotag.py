#! /usr/bin/env python
"""Convert inline barcodes into haplotag style ones"""
import os
import sys
import gzip
import argparse
from itertools import zip_longest

parser = argparse.ArgumentParser(
    prog = 'inline_to_haplotag.py',
    description = 'Moves inline linked read barcodes to read headers (OX:Z) and converts them into haplotag ACBD format (BX:Z). Barcodes must all be the same length.',
    usage = f"inline_to_haplotag.py -f <forward.fq.gz> -r <reverse.fq.gz> -b <barcodes.txt> -p <prefix>",
    exit_on_error = False
    )
parser.add_argument("-f", "--forward", required = True, type = str, help = "Forward reads of paired-end FASTQ file pair (gzipped)")
parser.add_argument("-r", "--reverse", required = True, type = str, help = "Reverse reads of paired-end FASTQ file pair (gzipped)")
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

def process_record(fw_entry, rv_entry, barcode_dict, bc_len):
    """convert the barcode to haplotag"""
    # [0] = header, [1] = seq, [2] = +, [3] = qual
    if fw_entry:
        bc_inline = fw_entry[1][:bc_len]
        bc_hap = barcode_dict.get(bc_inline, "A00C00B00D00")
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

bc_dict = {}
nucleotides = {'A','C','G','T'}
lengths = set()
# read in barcodes
opener = gzip.open if args.barcodes.lower().endswith('.gz') else open
mode = 'rt' if args.barcodes.lower().endswith('.gz') else 'r'
with opener(args.barcodes, mode) as bc_file:
    for line in bc_file:
        try:
            ATCG,ACBD = line.rstrip().split()
        except ValueError:
            sys.stderr.write(f"Invalid barcode entry: {line.rstrip()}\nExpected two entries: a nucleotide barcode and ACBD format barcode with a space/tab separatng them, e.g. ATATCAGA A01C22B13D93")
            sys.exit(1)
        if not set(ATCG).issubset(nucleotides):
            sys.stderr.write(f"Invalid barcode format: {ATCG}. Barcodes must be captial letters and only contain standard nucleotide values ATCG.\n")
            sys.exit(1)
        bc_dict[ATCG] = ACBD
        lengths.add(len(ATCG))
    if len(lengths) > 1:
        sys.stderr.write(f"Can only search sequences for barcodes of a single length, but multiple barcode legnths detected: " + ",".join([str(i) for i in lengths]))
    else:
        bc_len = lengths.pop()

# simultaneously iterate the forward and reverse fastq files
with gzip.open(args.forward, "r") as fw_i, gzip.open(args.reverse, "r") as rv_i,\
    gzip.open(f"{args.prefix}.R1.fq.gz", "wb", 6) as fw_out,\
    gzip.open(f"{args.prefix}.R2.fq.gz", "wb", 6) as rv_out:
    for fw_record, rv_record in zip_longest(iter_fastq_records(fw_i), iter_fastq_records(rv_i)),:
        new_fw, new_rv = process_record(fw_record, rv_record, bc_dict, bc_len)
        if new_fw:
            fw_out.write(new_fw.encode("utf-8"))
        if new_rv:
            rv_out.write(new_rv.encode("utf-8"))
