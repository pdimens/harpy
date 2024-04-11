import re
import sys
import gzip
import pysam
from itertools import chain

d = dict()
chromlast = False
alnfile = pysam.AlignmentFile(snakemake.input[0])
outfile = gzip.open(snakemake.output[0], "wb")
outfile.write(b"contig\tmolecule\treads\tstart\tend\tlength_inferred\tpercent_coverage\taligned_bp\tmindist\n")

def writestats(x,chr):
    """
    Writes molecule stats to a file. Only gets called when the chromosome changes.
    """
    for mi in x:
        x[mi]["inferred"] = x[mi]["end"] - x[mi]["start"] 
        x[mi]["mindist"] = max(0, x[mi]["mindist"])
        try:
            mi_aln = x[mi].get("alignments", 0)
            # converts start/end positions into ranges and merges ranges into a set, then gets the length of the unique set
            x[mi]["covered"] = round((len(set(chain(*[range(i,j+1) for i,j in mi_aln]))) - 1) / x[mi]["inferred"], 2)
        except:
            x[mi]["covered"] = 0
        outtext = f"{chr}\t{mi}\t" + "\t".join([str(x[mi][i]) for i in ["n", "start","end", "inferred", "covered", "bp", "mindist"]])
        outfile.write(outtext.encode() + b"\n")

for read in alnfile.fetch():
    chrm = read.reference_name
    bp   = read.query_alignment_length
    # check if the current chromosome is different from the previous one
    # if so, print the dict to file and empty it (a consideration for RAM usage)
    if chromlast != False and chrm != chromlast:
        writestats(d, chromlast)
        d = dict()
    
    chromlast = chrm
    
    if read.is_duplicate or read.is_unmapped:
        continue
    
    try:
        mi = read.get_tag("MI")
        bx = read.get_tag("BX")
        validBX = True
        # do a regex search to find X00 pattern in the BX
        if re.search("[ABCD]0{2,4}", bx):
            # if found, invalid
            bx = "invalidBX"
            mi = "invalidBX"
            validBX = False
    except:
        # There is no bx tag
        mi = "noBX"
        validBX = False
    aln = read.get_blocks()
    if not aln:
        # unaligned, skip
        continue

    # if invalid/absent BX, skip the distance stuff
    if mi in ["noBX", "invalidBX"]:
        if mi not in d.keys():
            d[mi] = {
                "start":  0,
                "end": 0,
                "bp":   bp,
                "n":    1,
                "lastpos" : 0,
                "mindist" : 0
            }
        else:
            d[mi]["bp"] += bp
            d[mi]["n"] += 1
        continue

    # logic to accommodate split reads 
    # start position of first alignment
    pos_start = aln[0][0]
    # end position of last alignment
    pos_end   = aln[-1][1]

    # create bx entry if it's not present
    if mi not in d.keys():
        d[mi] = {
            "start":  pos_start,
            "end": pos_end,
            "bp":   bp,
            "n":    1,
            "lastpos" : pos_end,
            "alignments" : aln,
            "mindist" : -1
        }
        continue

    # update the basic alignment info of the barcode
    d[mi]["alignments"] += aln
    d[mi]["bp"] += bp
    d[mi]["n"]  += 1
    # only if low < currentlow or high > currenthigh
    if pos_start < d[mi]["start"]:
        d[mi]["start"] = pos_start
    if pos_end > d[mi]["end"]:
        d[mi]["end"] = pos_end

    if read.is_reverse or (read.is_forward and not read.is_paired):
        # set the last position to be the end of current alignment
        d[mi]["lastpos"] = pos_end

    # distance from last alignment = current aln start - previous aln end
    dist = pos_start - d[mi]["lastpos"]

    # only calculate the minimum distance between alignments
    # if it's a forward read or an unpaired reverse read
    if read.is_forward or (read.is_reverse and not read.is_paired):
        if dist < d[mi]["mindist"] or d[mi]["mindist"] < 0:
            d[mi]["mindist"] = dist

outfile.close()
