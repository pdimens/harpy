#!/usr/bin/env python
f = open(snakemake.log[0], "w")
import sys
sys.stderr = sys.stdout = f

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
    for seq in bar_codes[code_letter]:
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

def get_read_codes(index_read, left_segment, right_segment):
    left  = index_read[0:6] # protocol-dependent
    right = index_read[7:]
    status = "found"
    if left in bar_codes[left_segment]:
        lc = bar_codes[left_segment][left]
    else:
        lc = get_min_dist(left, left_segment)
        status = "corrected"

    if right in bar_codes[right_segment]:
        rc = bar_codes[right_segment][right]
    else:
        rc = get_min_dist(right, right_segment)
        status = "corrected"

    if (lc == f"{left_segment}00" or rc == f"{right_segment}00"):
        status = "unclear"
    return rc, lc, status


qxrx = snakemake.params.qxrx
part = snakemake.params.part
schema = snakemake.input.schema
keep_unknown = snakemake.params.keep_unknown
outdir = snakemake.params.outdir
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
if keep_unknown:
    samples.append("_unknown_sample")
#create an array of files (one per sample) for writing
R1_output = {sample: open(f"{outdir}/{sample}.{part}.R1.fq", 'w') for sample in samples}
R2_output = {sample: open(f"{outdir}/{sample}.{part}.R2.fq", 'w') for sample in samples}
segments = {'A':'', 'B':'', 'C':'', 'D':''}
unclear_read_map={}
clear_read_map={}
with (
    pysam.FastxFile(r1) as R1,
    pysam.FastxFile(r2) as R2,
    pysam.FastxFile(i1, persist = False) as I1,
    pysam.FastxFile(i2, persist = False) as I2,
    open(snakemake.output.bx_info, 'w') as BC_log
):
    for r1_rec, r2_rec, i1_rec, i2_rec in zip(R1, R2, I1, I2):
        segments['A'], segments['C'], R1_status = get_read_codes(i1_rec.sequence, "C", "A")
        segments['B'], segments['D'], R2_status = get_read_codes(i2_rec.sequence, "D", "B")
        statuses = [R1_status, R2_status]
        BX_code = segments['A'] + segments['C'] + segments['B']+ segments['D']
        bc_tags = f"BX:Z:{BX_code}"
        if qxrx:
            bc_tags = f"RX:Z:{i1_rec.sequence}+{i2_rec.sequence}\tQX:Z:{i1_rec.quality}+{i2_rec.quality}\t{bc_tags}"
        r1_rec.comment += f"\t{bc_tags}"
        r2_rec.comment += f"\t{bc_tags}"
        # search sample name
        sample_name = samples_dict.get(segments[id_letter], "_unknown_sample")
        if sample_name == "_unknown_sample" and not keep_unknown:
            continue
        R1_output[sample_name].write(f"{r1_rec}\n")
        R2_output[sample_name].write(f"{r2_rec}\n")

        if "unclear" in statuses:
            if BX_code in unclear_read_map:
                unclear_read_map[BX_code] += 1
            else:
                unclear_read_map[BX_code] = 1
        else:        
            if "corrected" in statuses:
                if  BX_code in clear_read_map:
                    clear_read_map[BX_code][1] += 1
                else:
                    clear_read_map[BX_code] = [0,1] 
            else:
                if all(status == "found" for status in statuses):
                    if  BX_code in clear_read_map:
                        clear_read_map[BX_code][0] += 1
                    else:
                        clear_read_map[BX_code] = [1,0]          

    for sample_name in samples:
        R1_output[sample_name].close()
        R2_output[sample_name].close()

    BC_log.write("Barcode\tTotal_Reads\tCorrect_Reads\tCorrected_Reads\n")
    for code in clear_read_map:
        BC_log.write(f"{code}\t{sum(clear_read_map[code])}\t{clear_read_map[code][0]}\t{clear_read_map[code][1]}\n")
    for code in unclear_read_map:
        BC_log.write(f"{code}\t{unclear_read_map[code]}\t0\t0\n")
f.close()