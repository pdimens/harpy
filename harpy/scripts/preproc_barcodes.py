#! /usr/bin/env python3
import argparse
import json
import os
import pysam
import sys

def main():
    parser = argparse.ArgumentParser(
        prog='preproc-barcodes',
        description='[INTERNAL] script to rename and record demultiplexed linked read barcodes coming out of Pheniqs.',
        usage = "preproc-barcodes pheniqs.config.json sample.bam > out.bam"
        )
    parser.add_argument("json", type=str, help="Pheniqs configuration file that has the barcode definitions")
    parser.add_argument("bam", type=str, help="SAM/BAM file produced by Pheniqs")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    if not os.path.isfile(args.json):
        parser.error(f"{args.json} does not exist.")
    if not os.path.isfile(args.bam):
        parser.error(f"{args.bam} does not exist.")

    def reconstruct_barcode(nuc_bc: str):
        segments = nuc_bc.split("-")
        A = "A" + bc.get(segments[1], "00")
        B = "B" + bc.get(segments[2], "00")
        C = "C" + bc.get(segments[3], "00")
        D = "D" + stagger.get(segments[0], "00")
        valid = not any([A=="A00", B=="B00", C=="C00", D=="D00"])
        return f"{A}{C}{B}{D}", valid
        
    with open(args.json, 'r') as file:
        data = json.load(file)
    stagger = {}; bc = {}
    for k,v in data['decoder']['stagger']['codec'].items():
        stagger[v['barcode'][0]] = k.removeprefix("@")
    for k,v in data['decoder']['segment']['codec'].items():
        bc[v['barcode'][0]] = k.removeprefix("@")
    del data

    with (
        pysam.AlignmentFile(args.bam, require_index=False, check_sq = False) as bam,
        pysam.AlignmentFile(sys.stdout.buffer, 'w', template=bam, add_sam_header=False) as out
        ):
        for record in bam.fetch(until_eof=True):
            if record.has_tag("RX"):
                bx, vx = reconstruct_barcode(record.get_tag("RX"))
                record.set_tag("BX", bx, 'Z')
                record.set_tag("VX", vx, 'i')
            out.write(record)