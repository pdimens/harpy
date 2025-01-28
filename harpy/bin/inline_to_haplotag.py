#! /usr/bin/env python
"""Convert inline barcodes into haplotag style ones"""
import os
import sys
import gzip
import sqlite3
import argparse
from itertools import zip_longest

parser = argparse.ArgumentParser(
    prog = 'inline_to_haplotag.py',
    description = 'Moves inline linked read barcodes to read headers (OX:Z) and converts them into haplotag ACBD format (BX:Z). Barcodes must all be the same length.',
    usage = f"inline_to_haplotag.py -b <barcodes.txt> -p <prefix> FORWARD.fq.gz REVERSE.fq.gz",
    exit_on_error = False
    )
parser.add_argument("-p", "--prefix", required = True, type = str, help = "Prefix for outfile files (e.g. <prefix>.R1.fq.gz)")
parser.add_argument("-b", "--barcodes", required = True, type=str, help="Barcode conversion key file with format: ATCG<tab>ACBD")
parser.add_argument("forward", required = True, type = str, help = "Forward reads of paired-end FASTQ file pair (gzipped)")
parser.add_argument("reverse", required = True, type = str, help = "Reverse reads of paired-end FASTQ file pair (gzipped)")
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

def valid_record(fq_rec, FR):
    """fastq format sanity check"""
    if not (fq_rec[0].startswith("@") and fq_rec[2] == "+"):
        raise ValueError(f"Invalid FASTQ format for {FR} reads")

def insert_key_value(conn, key, value):
    """insert a key-value pair into sqlite database"""
    cursor = conn.cursor()
    cursor.execute('''
        INSERT OR REPLACE INTO kv_store (key, value) VALUES (?, ?)
    ''', (key, value))  # Use parameterized queries to avoid SQL injection
    conn.commit()

def get_value_by_key(conn, key):
    """retrieve a value by key from sqlite database"""
    cursor = conn.cursor()
    cursor.execute('''
        SELECT value FROM kv_store WHERE key = ?
    ''', (key,))
    result = cursor.fetchone()  # Fetch one row
    if result:
        # Return the value (first column of the result)
        return result[0]  
    else:
        # Return invalid ACBD haplotag if the key does not exist
        return "A00C00B00D00"  

def process_record(fw_entry, rv_entry, barcode_database, bc_len):
    """convert the barcode to haplotag"""
    # [0] = header, [1] = seq, [2] = +, [3] = qual
    if fw_entry:
        valid_record(fw_entry, "forward")
        bc_inline = fw_entry[1][:bc_len]
        bc_hap = get_value_by_key(barcode_database, bc_inline)
        fw_entry[0] = fw_entry[0].split()[0] + f"\tOX:Z:{bc_inline}\tBX:Z:{bc_hap}"
        fw_entry[1] = fw_entry[1][bc_len:]
        fw_entry[3] = fw_entry[3][bc_len:]
        _new_fw = "\n".join(fw_entry) + "\n"
        if rv_entry:
            valid_record(rv_entry, "reverse")
            rv_entry[0]  = rv_entry[0].split()[0] + f"\tOX:Z:{bc_inline}\tBX:Z:{bc_hap}"
            _new_rv = "\n".join(rv_entry) + "\n"
        else:
            _new_rv = None
    else:
        _new_fw = None
        # no forward read, therefor no barcode to search for
        if rv_entry:
            valid_record(rv_entry, "reverse")
            rv_entry[0] = rv_entry[0].split()[0] + "\tBX:Z:A00C00B00D00"
            _new_rv = "\n".join(rv_entry) + "\n"
        else:
            _new_rv = None
    return _new_fw, _new_rv

# Connect to an in-memory SQLite database
bc_db = sqlite3.connect(':memory:')
# Create the table to store key-value pairs
bc_db.cursor().execute('''
    CREATE TABLE kv_store (
        key TEXT PRIMARY KEY,
        value TEXT
    )
''')
bc_db.commit()

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
        
        insert_key_value(bc_db, ATCG, ACBD)
        lengths.add(len(ATCG))
    if len(lengths) > 1:
        sys.stderr.write("Can only search sequences for barcodes of a single length, but multiple barcode legnths detected: " + ",".join([str(i) for i in lengths]))
    else:
        bc_len = lengths.pop()

# simultaneously iterate the forward and reverse fastq files
with gzip.open(args.forward, "r") as fw_i, gzip.open(args.reverse, "r") as rv_i,\
    gzip.open(f"{args.prefix}.R1.fq.gz", "wb", 6) as fw_out,\
    gzip.open(f"{args.prefix}.R2.fq.gz", "wb", 6) as rv_out:
    record_F = []
    record_R = []
    for fw_record, rv_record in zip_longest(fw_i, rv_i):
        try:
            record_F.append(fw_record.decode().rstrip("\n"))
        except AttributeError:
            # if the file ends before the other one
            pass
        try:
            record_R.append(rv_record.decode().rstrip("\n"))
        except AttributeError:
            pass
        # sanity checks
        if len(record_F) == 4 or len(record_R) == 4:
            new_fw, new_rv = process_record(record_F, record_R, bc_db, bc_len)
            if new_fw:
                fw_out.write(new_fw.encode("utf-8"))
                record_F = []
            if new_rv:
                rv_out.write(new_rv.encode("utf-8"))
                record_R = []

bc_db.cursor().close()
