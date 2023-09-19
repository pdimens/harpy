#!/usr/bin/env python

import pysam
import sys
import re
import os.path

bam_in = sys.argv[1]

# regex for EXACTLY AXXCXXBXXDXX
haplotag = re.compile('^A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]$')
bam_pattern = re.compile("\.[bB][aA][mM]$", flags = re.IGNORECASE)

corename = re.sub(bam_pattern, "", os.path.basename(bam_in))

alnfile = pysam.AlignmentFile(bam_in)
if alnfile.header.get("RG")[0]['ID'] == corename:
    nameMismatch = 0
else:
    nameMismatch = 1

n_reads      = 0
noBX         = 0
badBX        = 0

for record in alnfile.fetch():
    n_reads += 1
    try:
        bx = record.get_tag("BX")
    except:
        # There is no bx tag
        noBX += 1
        continue
    # do a regex search to find AXXCXXBXXDXX pattern in the BX
    if not re.search(haplotag, bx):
        # malformed BX tag
        badBX += 1
        continue

alnfile.close()


values = [str(i) for i in [os.path.basename(bam_in), nameMismatch, n_reads, noBX, badBX]]
print("\t".join(values), file = sys.stdout) 