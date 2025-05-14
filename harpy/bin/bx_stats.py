#! /usr/bin/env python
"""calculate linked-read metrics from barcode information"""
import os
import re
import sys
import subprocess
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'bx_stats.py',
    description =
    """
    Calculates various linked-read molecule metrics from the input alignment file.
    The alignment file is expected to be in "standard" linked-read format, that is,
    the barcode is contained in the BX:Z tag and the barcode validation is stored
    in the VX:i tag as 0 (invalid) or 1 (valid). Metrics include (per molecule): 
    number of reads, position start, position end, length of molecule inferred from
    alignments, total aligned basepairs, total, length of inferred inserts, molecule
    coverage (%) based on aligned bases, molecule coverage (%) based on total inferred
    insert length. Input file **must be coordinate sorted**.
    """,
    usage = "bx_stats.py input.bam > output.gz",
    exit_on_error = False
    )

parser.add_argument('input', help = "Input coordinate-sorted bam/sam file.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
if not os.path.exists(args.input):
    parser.error(f"{args.input} was not found")

def writestats(x, writechrom, destination):
    """write to file the bx stats dictionary as a table"""
    for _mi in list(x.keys()):
        x[_mi]["inferred"] = x[_mi]["end"] - x[_mi]["start"]
        cov_bp = max(0, round(min(x[_mi]["bp"] / x[_mi]["inferred"], 1.0),4))
        x[_mi]["covered_bp"] = max(0, round(min(x[_mi]["bp"] / x[_mi]["inferred"], 1.0),4))
        x[_mi]["covered_inserts"] = max(0, round(min(x[_mi]["insert_len"] / x[_mi]["inferred"], 1.0), 4))
        outtext = f"{writechrom}\t{_mi}\t" + "\t".join([str(x[_mi][i]) for i in ["n", "start","end", "inferred", "bp", "insert_len", "covered_bp", "covered_inserts"]])
        destination.stdin.write(f"{outtext}\n".encode("utf-8"))
        # delete the entry after processing to ease up system memory
        del x[_mi]

with pysam.AlignmentFile(args.input) as alnfile:
    gzip = subprocess.Popen(["gzip"], stdin = subprocess.PIPE, stdout = sys.stdout)
    gzip.stdin.write(b"contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n")

    d = {}
    all_bx = set()
    LAST_CONTIG = None

    for read in alnfile.fetch(until_eof=True):
        chrom = read.reference_name
        # check if the current chromosome is different from the previous one
        # if so, print the dict to file and empty it (a consideration for RAM usage)
        if LAST_CONTIG and chrom != LAST_CONTIG:
            writestats(d, LAST_CONTIG, gzip)
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
        bp = read.query_alignment_length

        try:
            mi = read.get_tag("MI")
            bx = read.get_tag("BX")
            valid = bool(int(read.get_tag("VX")))
            # do a regex search to find X00 pattern in the BX
            if not valid:
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
        except KeyError:
            # There is no bx/MI tag
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

        # start position of first alignment
        ref_positions = [read.reference_start, read.reference_end]
        pos_start = min(ref_positions)
        # end position of last alignment
        pos_end   = max(ref_positions)
        if read.is_paired:
            if read.is_supplementary:
                # if it's a supplementary alignment, just use the alignment length
                isize = bp
            else:
                # by using max(), will either add 0 or positive TLEN to avoid double-counting
                isize = max(0, read.template_length)
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
            if read.is_forward:
                # +1 for a forward read, whether it is paired or not
                d[mi]["n"]  += 1
            elif read.is_reverse and not read.is_paired:
                # +1 for reverse only if it's unpaired, so the paired read doesn't count twice
                d[mi]["n"]  += 1
            d[mi]["bp"] += bp
            d[mi]["insert_len"] += isize
            d[mi]["start"] = min(pos_start, d[mi]["start"])
            d[mi]["end"] = max(pos_end, d[mi]["end"])

    # print the last entry
    writestats(d, LAST_CONTIG, gzip)
    # write comment on the last line with the total number of unique BX barcodes
    gzip.stdin.write(f"#total unique barcodes: {len(all_bx)}\n".encode("utf-8"))
