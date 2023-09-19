#!/usr/bin/env python

import pysam
import sys
import re
import os.path

fq_in = sys.argv[1]

bxz = re.compile('BX:Z:')
samspec = re.compile('[A-Z][A-Z]:[AifZHB]:')
haplotag = re.compile('A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]')

with pysam.FastxFile(fq_in) as fh:
    n_reads    = 0
    noBX       = 0
    badBX      = 0
    badSamSpec = 0
    for entry in fh:
        n_reads += 1
        # look for BX:Z: tag
        if bxz.match(entry.comment):
            # if AXXCXXBXXDXX format isnt found, it's a bad BX
            if not haplotag.search(entry.comment):
                badBX += 1
            for i in entry.comment.split():
                # if comments dont start with TAG:TYPE:, invalid SAM spec
                if not samspec.match(i):
                    badSamSpec += 1
        else:
            # missing BX:Z: tag
            noBX += 1

values = [str(i) for i in [os.path.basename(fq_in), n_reads, noBX, badBX, badSamSpec]]
print("\t".join(values), file = sys.stdout) 