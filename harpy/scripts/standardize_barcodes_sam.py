#! /usr/bin/env python3

import argparse
import sys
import re
import pysam

def bx_search(rec):
    if rec.has_tag('BX'):
        # presumably haplotagging
        if not rec.has_tag("VX"):
            rec.set_tag("VX", 0 if "00" in rec.get_tag("BX") else 1, value_type="i")
        return rec
    bx_pattern = r"(?:#\d+_\d+_\d+|\:[ATCGN]+)$"
    bx = re.search(bx_pattern, rec.query_name)
    if bx:
        barcode = bx[0]
        # set validation tag VX:i
        if barcode.startswith("#"):
            rec.set_tag("VX", 0 if "0" in barcode.split("_") else 1, value_type="i")
        else:
            rec.set_tag("VX", 0 if "N" in barcode else 1, value_type="i")
        rec.set_tag("BX", barcode[1:])
        rec.query_name = rec.query_name.removesuffix(barcode)
    return rec

def main():
    parser = argparse.ArgumentParser(
        prog='standardize_barcodes_sam',
        description='Convert barcode notation in input SAM file to BX:Z:BARCODE VX:i:0/1, where 0/1 is whether the barcode is valid (1) or invalid (0)',
        usage = "standardize_barcodes_sam input.sam > output.sam",
        )
    parser.add_argument('input_sam', nargs='?', type=argparse.FileType('r'), default= (None if sys.stdin.isatty() else sys.stdin), help = "Input SAM file")
    args = parser.parse_args()
    if not args.input_sam:
        parser.print_help(sys.stderr)
        sys.exit(1)

    try:
        with (
            pysam.AlignmentFile(args.input_sam, require_index=False) as SAM,
            pysam.AlignmentFile(sys.stdout.buffer, "w", template=SAM) as OUT
        ):
            for record in SAM.fetch(until_eof=True):
                OUT.write(bx_search(record))
    except (IOError, ValueError) as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)