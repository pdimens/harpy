#! /usr/bin/env python3

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
        if barcode.startswith("#"):
        # set validation tag VX:i
            rec.set_tag("VX", 0 if "0" in barcode.split("_") else 1, value_type="i")
        else:
            rec.set_tag("VX", 0 if "N" in barcode else 1, value_type="i")
        rec.set_tag("BX", barcode[1:])
        rec.query_name = rec.query_name.removesuffix(barcode)
    return rec

if not sys.stdin.isatty():
    # Input is being piped in via stdin
    input_data = sys.stdin
elif len(sys.argv) > 1:
    # First argument passed via command line
    input_data = sys.argv[1]

with pysam.AlignmentFile(input_data, require_index=False) as SAM, pysam.AlignmentFile(sys.stdout.buffer, "w", template=SAM) as OUT:
    for record in SAM.fetch(until_eof=True):
        OUT.write(bx_search(record))