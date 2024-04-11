#! /usr/bin/env python

## the ARGS ##
# 1 forward read
# 2 reverse read
# 3 barcode file
# 4 barcode conversion output file name
import gzip
import sys
from itertools import zip_longest, product

def process_record(fw_entry, rv_entry):
    # [0] = header, [1] = seq,[2] = +, [3] = qual
    bc10x = fw_entry[1][:16]
    bchap = bc_dict.get(bc10x, "A00C00B00D00")
    if not bchap:
        bchap = "".join(next(bc_generator))
        bc_dict[bc10x] = bchap
    new_fw  = fw_entry[0].split()[0] + f"\tTX:Z:{bc10x}\tBX:Z:{bchap}\n"
    new_fw += fw_entry[1][16:] + "\n"
    new_fw += fw_entry[2] + "\n"
    new_fw += fw_entry[3][16:] + "\n"
    new_rv  = rv_entry[0].split()[0] + f"\tTX:Z:{bc10x}\tBX:Z:{bchap}\n"
    new_rv += "\n".join(rv_entry[1:3])
    return new_fw, new_rv

bc_range = [f"{i}".zfill(2) for i in range(1,97)]
bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)

bc_dict = dict()

# read in barcodes
with open(snakemake.input[2], "r") as bc_file:
    while True:
        # Get next line from file
        if snakemake.input[2].endswith("gz"):
            line = bc_file.readline().decode()
        else:
            line = bc_file.readline()
        # if line is empty
        # end of file is reached
        if not line:
            break
        bc = line.rstrip("\n").split()
        _10x = str(bc[0])
        bc_dict[_10x] = None

# simultaneously iterate the forward and reverse fastq files
fw_reads = snakemake.input[0]
rv_reads = snakemake.input[1]

fw_out = gzip.open(snakemake.output[0], "w")
rv_out = gzip.open(snakemake.output[1], "w")

with gzip.open(fw_reads, "r") as fw_i, gzip.open(rv_reads, "r") as rv_i:
    # store the fq records here
    fw_record = []
    rv_record = []
    i = 0
    for fw, rv in zip_longest(fw_i, rv_i):
        fw_line = fw.decode().rstrip("\n")
        rv_line = rv.decode().rstrip("\n")
        if fw_line.startswith("@") and i > 0:
            itera += 1
            # process the full record
            new_fw, new_rv = process_record(fw_record, rv_record)
            # write new record to files
            fw_out.write(new_fw.encode())
            rv_out.write(new_rv.encode())
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

with open(snakemake.log[0], "w") as bc_out:
    for i in bc_dict:
        bc_out.write(i + "\t" + bc_dict[i] + "\n") if bc_dict[i] else None