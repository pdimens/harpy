#! /usr/bin/env python3

# THIS SCRIPT IS INCOMPLETE AND CURRENTLY REDUNDANT WITH writeBED.pl. Might replace it one day.
file = "~/ShadHap1_C01.molsize"

file = """
CM031716.1      10045943        10045962       19       A80C01B63D46    1
CM031716.1      1009539 1009558 19      A10C01B13D71    1
CM031716.1      1013338 1013357 19      A07C01B43D04    1
CM031716.1      1015177 1015196 19      A83C01B88D63    1
CM031716.1      1019277 1019296 19      A80C01B20D91    1
CM031716.1      102313  102332  19      A02C01B50D53    1
CM031716.1      10234632        10234651       19       A37C01B83D04    1
CM031716.1      10249773        10249792       19       A48C01B18D92    1
CM031716.1      1030023 1030042 19      A84C01B35D33    1
CM031716.1      1035577 1035596 19      A38C01B75D71    1
"""

d = dict()

# "samtools view $ARGV[0] -F 1024

with open(file, "r") as f:
    while True:
        line = f.readline().split()
        if not line:
            break
        chrm = line[0]
        bx = line[4]
        lw = int(line[1])
        hi = int(line[2])
        # create chromosome key if not present
        # populate it with the bx stats
        if chrm not in d.keys():
            d[chrm] = {
                bx : {
                "low":  lw,
                "high": hi
                }    
            }
        # create bx stats if it's not present
        elif bx not in d[chrm].keys():
            d[chrm] = {
                bx : {
                "low":  lw,
                "high": hi
                }    
            }
        # if BX is present for this chrm, update
        # only if low < currentlow or high > currenthigh
        else:
            if lw < d[chrm][bx]["low"]:
                d[chrm][bx]["low"] = lw
            if hi > d[chrm][bx]["high"]:
                d[chrm][bx]["high"] = hi