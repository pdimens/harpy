#!/usr/bin/env python

import os
import sys
import gzip
import argparse
import subprocess
import pysam
from Levenshtein import distance

parser = argparse.ArgumentParser(
    prog = 'dmux_gen1.py',
    description =
    """
    Demultiplex raw linked-reads that were created using the original Meier et al. haplotagging protocol.
    All inputs are required. The schema should have rows that look like "Sample_17<tab>C17". The --bx-{acbd}
    inputs should have rows that look like "A01<tab>ACGGAA" or "D22<tab>ACGCC" (without quotes), where the first column has
    the segment corresponding to the letter of the argument ("A" for --bx-a, "C" for --bx-c, etc.).
    """,
    usage = "demux_gen1.py --r1 data.R1.fq.gz --r2 data.R2.fq.gz --i1 data.I1.fq.gz --i2 data.I2.fq.gz --schema samples.schema --bx-{acbd} ...",
    exit_on_error = False
    )

parser.add_argument('--r1', required = True, type = str, help = "Forward reads, in fastq format")
parser.add_argument('--r2', required = True, type = str, help = "Reverse reads, in fastq format")
parser.add_argument('--i1', required = True, type = str, help = "Forward index reads, in fastq format")
parser.add_argument('--i2', required = True, type = str, help = "Reverse index reads, in fastq format")
parser.add_argument('--schema', required = True, type = str, help = "Sample demultiplexing schema of sample<tab>id_segment")
parser.add_argument('--bx-a', type = str, required = True, help = "file of nucleotide_barcode<tab>A_segement")
parser.add_argument('--bx-b', type = str, required = True, help = "file of nucleotide_barcode<tab>B_segement")
parser.add_argument('--bx-c', type = str, required = True, help = "file of nucleotide_barcode<tab>C_segement")
parser.add_argument('--bx-d', type = str, required = True, help = "file of nucleotide_barcode<tab>D_segement")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
err = []
for arg_name, arg_val in vars(args).items():
    if not os.path.exists(arg_val):
        err.append(arg_val)
if err:
    parser.error("Input files were not found:\n" + ", ".join(err))

def read_barcodes(file_path, segment):
    """Read and parse input barcode (segment) file of segment<tab>sequence"""
    data_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            try:
                code, seq = line.rstrip().split()
                if code[0].upper() != segment:
                    parser.error(f"Segments in {file_path} are expected to begin with {segment}, but begin with {code[0].upper()}")
                data_dict[seq] = code
            except ValueError:
                # skip rows without two columns
                continue
    return data_dict

def read_schema(file_path):
    """Read and parse schema file of sample<tab>id_segment"""
    # one sample can have more than one code
    # {segment : sample}
    data_dict = {}
    code_letters = set() #codes can be Axx, Bxx, Cxx, Dxx
    with open(file_path, 'r') as file:
        for line in file:
            try:
                sample, segment_id = line.rstrip().split()
                data_dict[segment_id] = sample
                code_letters.add(segment_id[0])
            except ValueError:
                # skip rows without two columns
                continue
    if not code_letters:
        print(f"Schema file {os.path.basename(file_path)} has no valid rows for demultiplexing. Rows should be sample<tab>segment, e.g. sample_01<tab>C75")
        sys.exit(1)
    if len(code_letters) > 1:
        print(f"Schema file {os.path.basename(file_path)} has sample IDs in more than one linked-read segment. All sample IDs should be in a single segment, such as C or D")
        sys.exit(1)
    id_letter = code_letters.pop()
    return id_letter, data_dict

def get_min_dist(needle, code_letter):
    minDist = 999
    nbFound = 0
    minSeq =""
    for seq in bar_codes[code_letter].keys():
        d = distance(needle, seq)
        if (d < minDist):
            minDist = d
            nbFound = 1
            minSeq= seq
        elif (d == minDist):
            nbFound += 1   
    if (nbFound>1):
        code_min_dist = f"{code_letter}00"
    else:
        code_min_dist =  bar_codes[code_letter][minSeq]   
    return code_min_dist

def get_read_codes(I, codeL, codeR):
  left  = I[0:6] # change depending on protocol
  right = I[7:]
  status = "found"
  if left in bar_codes[codeL]:
    lc = bar_codes[codeL][left]
  else:
    lc = get_min_dist(left, codeL)
    status = "corrected"

  if right in bar_codes[codeR]:
    rc = bar_codes[codeR][right]
  else:
    rc = get_min_dist(right, codeR)
    status = "corrected"

  if (lc == f"{codeL}00" or rc == f"{codeR}00"):
    status = "unclear"
  return rc, lc, status

bar_codes = {
    "A" : read_barcodes(args.bx_a, "A"),
    "B" : read_barcodes(args.bx_b, "B"),
    "C" : read_barcodes(args.bx_c, "C"),
    "D" : read_barcodes(args.bx_d, "D"),
}

#read schema
id_letter, samples_dict = read_schema(args.schema)
samples = list(set(samples_dict.values()))
samples.append("unknown_data")
#create an array of files (one per sample) for writing
R1_output = {sample: gzip.open(sample + ".R1.fastq.gz", 'wb', compresslevel = 6) for sample in samples}
R2_output = {sample: gzip.open(sample + ".R2.fastq.gz", 'wb', compresslevel = 6) for sample in samples}

segments = {'A':'', 'B':'', 'C':'', 'D':''}
unclear_read_map={}
clear_read_map={}
with (
    pysam.FastxFile(args.r1) as R1,
    pysam.FastxFile(args.r2) as R2,
    pysam.FastxFile(args.i1, persist = False) as I1,
    pysam.FastxFile(args.i2, persist = False) as I2,
    open('clearBC.log', 'w') as clearBC_log,
    open('unclearBC.log', 'w') as unclearBC_log
):
    for r1_rec, r2_rec, i1_rec, i2_rec in zip(R1, R2, I1, I2):
        segments['A'], segments['C'], statusR1 = get_read_codes(i1_rec.sequence, "C", "A")
        segments['B'], segments['D'], statusR2 = get_read_codes(i2_rec.sequence, "D", "B")
        BX_code = segments['A'] + segments['C'] + segments['B']+ segments['D']
        # search sample name
        sample_name = samples_dict.get(segments[id_letter], "unknown_data")
        bc_tags = [
            f"RX:Z:{i1_rec.sequence}+{i2_rec.sequence}",
            f"QX:Z:{i1_rec.quality}+{i2_rec.quality}",
            f"BX:Z:{BX_code}"
        ]
        r1_rec.comment += "\t" + "\t".join(bc_tags)
        r2_rec.comment += "\t" + "\t".join(bc_tags)
        R1_output[sample_name].write(f"{r1_rec}\n".encode("utf-8"))
        R2_output[sample_name].write(f"{r2_rec}\n".encode("utf-8"))

        if (statusR1 == "unclear" or statusR2 == "unclear"):
            if BX_code in unclear_read_map:
                unclear_read_map[BX_code] += 1
            else:
                unclear_read_map[BX_code] = 1
        else:        
            if (statusR1 == "corrected" or statusR2 == "corrected"):
                if  BX_code in clear_read_map:
                    clear_read_map[BX_code][1] += 1
                else:
                    clear_read_map[BX_code] = [0,1] 
            else:
                if (statusR1 == "found" and statusR2 == "found"): 
                    if  BX_code in clear_read_map:
                        clear_read_map[BX_code][0] += 1
                    else:
                        clear_read_map[BX_code] = [1,0]            

    for sample_name in samples:
        R1_output[sample_name].close()
        R2_output[sample_name].close()

    clearBC_log.write("Barcode\tCorrect reads\tCorrected reads\n" )
    for code in clear_read_map:
        clearBC_log.write(code+"\t"+"\t".join(str(x) for x in clear_read_map[code])+"\n") 

    unclearBC_log.write("Barcode\tReads\n")
    for code in unclear_read_map:
        unclearBC_log.write(code +"\t"+str(unclear_read_map [code])+"\n")
