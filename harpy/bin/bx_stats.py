#! /usr/bin/env python
"""calculate linked-read metrics from barcode information"""
import re
import sys
import gzip
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'bx_stats.py',
    description =
    """
    Calculates various linked-read molecule metrics from the input alignment file.
    Metrics include (per molecule): number of reads, position start, position end,
    length of molecule inferred from alignments, total aligned basepairs, total,
    length of inferred inserts, molecule coverage (%) based on aligned bases, molecule
    coverage (%) based on total inferred insert length.
    Input file MUST BE COORDINATE SORTED.
    """,
    usage = "bx_stats.py -o output.gz input.bam",
    exit_on_error = False
    )

parser.add_argument('-o', '--output', help = "Gzipped tab-delimited file of metrics.")
parser.add_argument('input', help = "Input coordinate-sorted bam/sam file. If bam, a matching index file should be in the same directory.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


alnfile = pysam.AlignmentFile(args.input)
outfile = gzip.open(args.output, "wb", 6)
outfile.write(b"contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n")

d = {}
all_bx = set()
LAST_CONTIG = None

def writestats(x, writechrom):
    """write to file the bx stats dictionary as a table"""
    for _mi in x:
        x[_mi]["inferred"] = x[_mi]["end"] - x[_mi]["start"]
        try:
            x[_mi]["covered_bp"] = round(min(x[_mi]["bp"] / x[_mi]["inferred"], 1.0),4)
        except:
            x[_mi]["covered_bp"] = 0
        try:
            x[_mi]["covered_inserts"] = round(min(x[_mi]["insert_len"] / x[_mi]["inferred"], 1.0), 4)
        except:
            x[_mi]["covered_inferred"] = 0
        outtext = f"{writechrom}\t{_mi}\t" + "\t".join([str(x[_mi][i]) for i in ["n", "start","end", "inferred", "bp", "insert_len", "covered_bp", "covered_inserts"]])
        outfile.write(outtext.encode() + b"\n")

for read in alnfile.fetch():
    chrom = read.reference_name
    # check if the current chromosome is different from the previous one
    # if so, print the dict to file and empty it (a consideration for RAM usage)
    if LAST_CONTIG and chrom != LAST_CONTIG:
        writestats(d, LAST_CONTIG)
        d = {}
    LAST_CONTIG = chrom
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
        bx = read.get_tag("BX")
        # do a regex search to find X00 pattern in the BX
        if re.search("[ABCD]0{2,4}", bx):
            # if found, invalid
            if "invalidBX" not in d:
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
        # add valid bx to set of all unique barcodes
        # remove the deconvolve hyphen, if present
        all_bx.add(bx.split("-")[0])
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
        #isize = min(bp, abs(read.template_length))
        # only count the bp of the first read of paired end reads
        #bp = bp if read.is_read1 else 0
        bp = min(abs(read.template_length), read.infer_query_length())
    elif read.is_supplementary:
        # if it's a supplementary alignment, just use the alignment length
        isize = bp
    else:
        # if it's unpaired, use the TLEN or query length, whichever is bigger
        isize = max(abs(read.template_length), read.infer_query_length())
    # create bx entry if it's not present
    if mi not in d:
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
writestats(d, LAST_CONTIG)
# write comment on the last line with the total number of unique BX barcodes
outfile.write(f"#total unique barcodes: {len(all_bx)}\n".encode())
outfile.close()
