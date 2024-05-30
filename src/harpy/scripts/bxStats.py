import re
import sys
import gzip
import pysam

alnfile = pysam.AlignmentFile(snakemake.input[0])
outfile = gzip.open(snakemake.output[0], "wb", 6)
outfile.write(b"contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n")

d = dict()
chromlast = None

def writestats(x, contig):
    for mi in x:
        x[mi]["inferred"] = x[mi]["end"] - x[mi]["start"] 
        try:
            x[mi]["covered_bp"] = min(x[mi]["bp"] / x[mi]["inferred"], 1.0)
        except:
            x[mi]["covered_bp"] = 0
        try:
            x[mi]["covered_inserts"] = min(x[mi]["insert_len"] / x[mi]["inferred"], 1.0)
        except:
            x[mi]["covered_inferred"] = 0
        outtext = f"{chromlast}\t{mi}\t" + "\t".join([str(x[mi][i]) for i in ["n", "start","end", "inferred", "bp", "insert_len", "covered_bp", "covered_inserts"]])
        outfile.write(outtext.encode() + b"\n")

for read in alnfile.fetch():
    chrom = read.reference_name
    # check if the current chromosome is different from the previous one
    # if so, print the dict to file and empty it (a consideration for RAM usage)
    if chromlast and chrom != chromlast:
        writestats(d, chromlast)
        d = dict()
    chromlast = chrom
    # skip duplicates, unmapped, and secondary alignments
    if read.is_duplicate or read.is_unmapped or read.is_secondary:
        continue
    # skip chimeric alignments that map to different contigs
    if read.is_supplementary and read.reference_name != read.next_reference_name:
        continue
    if not read.get_blocks():
        continue

    # numer of bases aligned
    bp = read.reference_length

    try:
        mi = read.get_tag("MI")
        # do a regex search to find X00 pattern in the BX
        if re.search("[ABCD]0{2,4}", read.get_tag("BX")):
            # if found, invalid
            if "invalidBX" not in d.keys():
                d["invalidBX"] = {
                    "start":  0,
                    "end": 0,
                    "bp":   bp,
                    "insert_len" : 0,
                    "n":    1,
                }
            else:
                d["invalidBX"]["bp"] += bp
                d["invalidBX"]["n"] += 1
            continue
    except:
        # There is no bx/MI tag
        continue

    # start position of first alignment
    ref_positions = [read.reference_start, read.reference_end]
    pos_start = min(ref_positions)
    # end position of last alignment
    pos_end   = max(ref_positions)
    if read.is_paired:
        # by using max(), will either add 0 or positive TLEN to avoid double-counting
        isize = max(0, read.template_length)
        # only count the bp of the first read of paired end reads
        bp = bp if read.is_read1 else 0
    elif read.is_supplementary:
        # if it's a supplementary alignment, just use the alignment length
        isize = bp
    else:
        # if it's unpaired, use the TLEN or query length, whichever is bigger
        isize = max(abs(read.template_length), read.infer_query_length())
    # create bx entry if it's not present
    if mi not in d.keys():
        d[mi] = {
            "start": pos_start,
            "end": pos_end,
            "bp": bp,
            "insert_len": isize,
            "n": 1,
        }
    else:
        # update the basic alignment info of the molecule
        d[mi]["bp"] += bp
        d[mi]["n"]  += 1
        d[mi]["insert_len"] += isize
        d[mi]["start"] = min(pos_start, d[mi]["start"])
        d[mi]["end"] = max(pos_end, d[mi]["end"])

# print the last entry
writestats(d, chrom)
outfile.close()