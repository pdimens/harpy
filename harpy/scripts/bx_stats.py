#! /usr/bin/env python
"""calculate linked-read metrics from barcode information"""
import argparse
import os
import numpy as np
import pysam
import sys

class ReadCloud():
    '''
    A class to store the relevant information of alignment records that have the same `BX` barcode tag.
    '''
    def __init__(self, valid: bool = True):
        self.positions: list[list[int]] = []
        self.bp: list[int] = []
        self.inserts : list[int] = []
        self.count : list[bool] = []
        self.valid = valid
        self.barcode = ""

    def add(self, record: pysam.AlignedSegment):
        '''add a pysam alignment record to the read cloud, keeping only the relevant info'''
        if self.valid:
            self.positions.append([record.reference_start, record.reference_end])
            self.inserts.append(insert_size(record))
            self.barcode = record.get_tag("BX")
        self.bp.append(record.query_alignment_length)
        self.count.append(record.is_read1 or not record.is_paired)

    def deconvolve(self, chrom, cutoff) -> str:
        '''
        Process the read cloud and deconvolute using `cutoff` (if it's >0). Deconvolution
        appends `-N` to the barcodes where `N` is an integer (e.g. `-1`, `-2`). Returns
        a string of barcode row(s) and the associated stats.
        '''
        if not self.valid:
            return "{chrom}\t-1\t" + self.stats(0, 0, 0, sum(self.bp), sum(self.count))
        result = ""
        # sort alignment extrema by leftmost position and sort the subsequent info the same way
        sort_values = [sublist[0] for sublist in self.positions]
        sorted_indices = np.argsort(sort_values)
        self.positions = [self.positions[i] for i in sorted_indices]
        self.bp = [self.bp[i] for i in sorted_indices]
        self.inserts = [self.inserts[i] for i in sorted_indices]
        self.count = [self.count[i] for i in sorted_indices]

        # instantiate with the first value
        deconv = 0
        start = self.positions[0][0]
        end = self.positions[0][1]
        insert = self.inserts[0]
        bp = self.bp[0]
        count = int(self.count[0])

        for idx in range(1, len(self.positions)):
            prev_end = self.positions[idx - 1][1]
            curr_start = self.positions[idx][0]
            curr_end = self.positions[idx][1]

            # check if this alignment should be merged with current molecule
            if cutoff == 0 or (curr_start - prev_end) <= cutoff:
                # merge into current molecule
                end = max(end, curr_end)
                insert += self.inserts[idx]
                bp += self.bp[idx]
                count += self.count[idx]
            else:
                # gap exceeds cutoff - write current molecule and start new one
                BC = f"{self.barcode}-{deconv}" if deconv > 0 else self.barcode
                result += f"{chrom}\t{BC}\t" + self.stats(start, end, insert, bp, count)

                # start new molecule with current alignment
                deconv += 1
                start = self.positions[idx][0]
                end = curr_end
                insert = self.inserts[idx]
                bp = self.bp[idx]
                count = int(self.count[idx])

        # write final molecule
        BC = f"{self.barcode}-{deconv}" if deconv > 0 else self.barcode
        result += f"{chrom}\t{BC}\t" + self.stats(start, end, insert, bp, count)
        return result

    def stats(self, start, end, insert, bp, count) -> str:
        '''
        Calculate coverage and return a string of count, start position,
        end position, inferred length, bases aligned, insert length,
        coverage based on aligned bases, and coverage based on inferred inserts.
        '''
        #columns = reads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n"
        inferred = end - start
        try:
            cov_bp = max(0, round(min(bp / inferred, 1.0),5))
            cov_ins = max(0, round(min( insert / inferred, 1.0), 5))
        except ZeroDivisionError:
            cov_bp = 0
            cov_ins = 0
        return f"{count}\t{start}\t{end}\t{inferred}\t{bp}\t{insert}\t{cov_bp}\t{cov_ins}\n"

def writestats(x: dict, thresh, writechrom):
    """write to file the bx stats dictionary as a table"""
    for _mi in list(x.keys()):
        sys.stdout.write(x[_mi].deconvolve(writechrom, thresh))
        # delete the entry after processing to ease up system memory
        del x[_mi]

def insert_size(rec) -> int:
    '''Calculate the insert size'''
    # start position of first alignment
    if rec.is_paired:
        if rec.is_supplementary:
            # if it's a supplementary alignment, just use the alignment length
            isize = rec.query_alignment_length
        else:
            # by using max(), will either add 0 or positive TLEN to avoid double-counting
            isize = max(0, rec.template_length)
    else:
        # if it's unpaired, use the TLEN or query length, whichever is bigger
        isize = max(abs(rec.template_length), rec.infer_query_length())
    return isize

def main():
    parser = argparse.ArgumentParser(
        prog = 'bx_stats',
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
        usage = "bx-stats <-d DISTANCE> input.bam > output.gz",
        exit_on_error = False
        )

    parser.add_argument('-d', '--distance-threshold', type = int, default=0, help = "Calculate statistics assuming this distance threshold for linking alignments sharing a barcode")
    parser.add_argument('input', help = "Input coordinate-sorted bam/sam file.")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not os.path.exists(args.input):
        parser.error(f"{args.input} was not found")
    dist_thresh = max(0, args.distance_threshold)
    sys.stdout.write("contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n")
    with pysam.AlignmentFile(args.input) as alnfile:
        d = {}
        LAST_CONTIG = None

        for read in alnfile.fetch(until_eof=True):
            chrom = read.reference_name
            # check if the current chromosome is different from the previous one
            # if so, print the dict to file and empty it (a consideration for RAM usage)
            if LAST_CONTIG and chrom != LAST_CONTIG:
                writestats(d, dist_thresh, LAST_CONTIG)
                d = {}
            LAST_CONTIG = chrom
            # skip duplicates, unmapped, and secondary alignments
            if read.is_duplicate or read.is_unmapped or read.is_secondary:
                continue
            # skip chimeric alignments that map to different contigs
            if read.is_supplementary and read.reference_name != read.next_reference_name:
                continue
            if not read.get_blocks():
                # unaligned, skip it
                LAST_CONTIG = chrom
                continue

            try:
                bx = read.get_tag("BX")
                if read.get_tag("VX") == 0:
                    # VX:i:0 is invalid
                    raise KeyError
            except KeyError:
                # There is either no BX or no/invalid VX tag
                if "invalid" not in d:
                    d["invalid"] = ReadCloud(valid = False) 
                d["invalid"].add(read)
                LAST_CONTIG = chrom
                continue

            if bx not in d:
                d[bx] = ReadCloud()

            d[bx].add(read)
            LAST_CONTIG = chrom

        # print the last entry
        writestats(d, dist_thresh, LAST_CONTIG)
