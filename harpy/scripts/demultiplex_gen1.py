#!/usr/bin/env python

import os
import sys
import pysam
from Levenshtein import distance

def read_barcodes(file_path, segment):
    """Read and parse input barcode (segment) file of segment<tab>sequence"""
    data_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            try:
                code, seq = line.rstrip().split()
                if code[0].upper() != segment:
                    sys.stderr.write(f"Segments in {file_path} are expected to begin with {segment}, but begin with {code[0].upper()}\n")
                    sys.exit(1)
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
    # codes can be Axx, Bxx, Cxx, Dxx
    code_letters = set()
    with open(file_path, 'r') as file:
        for line in file:
            try:
                sample, segment_id = line.rstrip().split()
                data_dict[segment_id] = sample
                code_letters.add(segment_id[0])
            except ValueError:
                # skip rows without two columns
                continue
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
  left  = I[0:6] # protocol-dependant
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

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    outdir = snakemake.params.outdir
    qxrx = snakemake.params.qxrx
    schema = snakemake.input.schema
    r1 = snakemake.input.R1
    r2 = snakemake.input.R2
    i1 = snakemake.input.I1
    i2 = snakemake.input.I2
    bx_a = snakemake.input.segment_a
    bx_b = snakemake.input.segment_b
    bx_c = snakemake.input.segment_c
    bx_d = snakemake.input.segment_d
    bar_codes = {
        "A" : read_barcodes(bx_a, "A"),
        "B" : read_barcodes(bx_b, "B"),
        "C" : read_barcodes(bx_c, "C"),
        "D" : read_barcodes(bx_d, "D"),
    }

    #read schema
    id_letter, samples_dict = read_schema(schema)
    samples = list(set(samples_dict.values()))
    samples.append("unidentified_data")
    #create an array of files (one per sample) for writing
    R1_output = {sample: open(f"{outdir}/{sample}.R1.fq", 'w') for sample in samples}
    R2_output = {sample: open(f"{outdir}/{sample}.R2.fq", 'w') for sample in samples}

    segments = {'A':'', 'B':'', 'C':'', 'D':''}
    unclear_read_map={}
    clear_read_map={}
    with (
        pysam.FastxFile(r1) as R1,
        pysam.FastxFile(r2) as R2,
        pysam.FastxFile(i1, persist = False) as I1,
        pysam.FastxFile(i2, persist = False) as I2,
        open(snakemake.output.valid, 'w') as clearBC_log,
        open(snakemake.output.invalid, 'w') as unclearBC_log
    ):
        for r1_rec, r2_rec, i1_rec, i2_rec in zip(R1, R2, I1, I2):
            segments['A'], segments['C'], statusR1 = get_read_codes(i1_rec.sequence, "C", "A")
            segments['B'], segments['D'], statusR2 = get_read_codes(i2_rec.sequence, "D", "B")
            BX_code = segments['A'] + segments['C'] + segments['B']+ segments['D']
            bc_tags = f"BX:Z:{BX_code}"
            if qxrx:
                bc_tags = f"RX:Z:{i1_rec.sequence}+{i2_rec.sequence}\tQX:Z:{i1_rec.quality}+{i2_rec.quality}\t{bc_tags}"
            r1_rec.comment += f"\t{bc_tags}"
            r2_rec.comment += f"\t{bc_tags}"
            # search sample name
            sample_name = samples_dict.get(segments[id_letter], "unidentified_data")
            R1_output[sample_name].write(f"{r1_rec}\n")
            R2_output[sample_name].write(f"{r2_rec}\n")

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
