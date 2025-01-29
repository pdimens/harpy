#! /usr/bin/env python
"""Convert inline barcodes into haplotag style ones. Assumes reads are properly paired."""
import os
import sys
import gzip
import sqlite3
import argparse
from itertools import zip_longest
import pysam

parser = argparse.ArgumentParser(
    prog = 'inline_to_haplotag.py',
    description = 'Moves inline linked read barcodes to read headers (OX:Z) and converts them into haplotag ACBD format (BX:Z). Barcodes must all be the same length.',
    usage = f"inline_to_haplotag.py -b <barcodes.txt> -p <prefix> FORWARD.fq.gz REVERSE.fq.gz",
    exit_on_error = False
    )
parser.add_argument("-p", "--prefix", required = True, type = str, help = "Prefix for outfile files (e.g. <prefix>.R1.fq.gz)")
parser.add_argument("-b", "--barcodes", required = True, type=str, help="Barcode conversion key file with format: ATCG<tab>ACBD")
parser.add_argument("forward", type = str, help = "Forward reads of paired-end FASTQ file pair (gzipped)")
parser.add_argument("reverse", type = str, help = "Reverse reads of paired-end FASTQ file pair (gzipped)")
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

def store_bc_as_sql(conn, bclist):
    """insert a key-value pair into sqlite database"""
    cursor = conn.cursor()
    cursor.executemany('INSERT INTO kv_store (key, value) VALUES (?, ?)', bclist)
#    cursor.execute('''
#        INSERT OR REPLACE INTO kv_store (key, value) VALUES (?, ?)
#    ''', (key, value))  # Use parameterized queries to avoid SQL injection
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

def process_record(fw_rec, rv_rec, barcode_database, bc_len):
    """convert the barcode to haplotag"""
    if fw_rec and len(fw_rec.sequence) > bc_len:
        bc_inline = fw_rec.sequence[:bc_len]
        bc_hap = get_value_by_key(barcode_database, bc_inline)
        fw_rec.comment = fw_rec.comment.split()[0] + f"\tOX:Z:{bc_inline}\tBX:Z:{bc_hap}"
        fw_rec.sequence = fw_rec.sequence[bc_len:]
        fw_rec.quality = fw_rec.quality[bc_len:]
        _new_fw = str(fw_rec) + "\n"
        if rv_rec:
            rv_rec.comment = rv_rec.comment.split()[0] + f"\tOX:Z:{bc_inline}\tBX:Z:{bc_hap}"
            _new_rv = str(rv_rec) + "\n"
        else:
            _new_rv = None
    else:
        _new_fw = None
        # no forward read, therefore no barcode to search for
        if rv_rec:
            rv_rec.comment = rv_rec.comment.split()[0] + "\tBX:Z:A00C00B00D00"
            _new_rv = str(rv_rec) + "\n"
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
bc_len = None

# read in barcodes
opener = gzip.open if args.barcodes.lower().endswith('.gz') else open
mode = 'rt' if args.barcodes.lower().endswith('.gz') else 'r'
bc_mem_data = []
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
        bc_mem_data.append((ATCG, ACBD))
        if bc_len and len(ATCG) != bc_len:
            sys.stderr.write("Can only search sequences for barcodes of a single length, but multiple barcode legnths detected.")
            sys.exit(1)
        else:
            bc_len = len(ATCG)

store_bc_as_sql(bc_db, bc_mem_data)
# flush this list out of memory b/c it's already stored in the sql table
del bc_mem_data

# simultaneously iterate the forward and reverse fastq files
with (
    pysam.FastxFile(args.forward) as fw,
    pysam.FastxFile(args.reverse) as rv,
    gzip.open(f"{args.prefix}.R1.fq.gz", "wb", 6) as fw_out,
    gzip.open(f"{args.prefix}.R2.fq.gz", "wb", 6) as rv_out,
):
    for fw_record, rv_record in zip_longest(fw, rv):
        new_fw, new_rv = process_record(fw_record, rv_record, bc_db, bc_len)
        if new_fw:
            fw_out.write(new_fw.encode("utf-8"))
        if new_rv:
            rv_out.write(new_rv.encode("utf-8"))

bc_db.cursor().close()
