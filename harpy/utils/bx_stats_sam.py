import click
import numpy as np
import pysam
import sys

class ReadCloud():
    '''
    A class to store the relevant information of alignment records that have the same `BX` barcode tag.
    '''
    def __init__(self, valid: bool = True):
        self.chromosome: str = ""
        self.positions: list[list[int]] = []
        self.bp: list[int] = []
        self.inserts : list[int] = []
        self.count : list[bool] = []
        self.valid: bool = valid
        self.barcode: str = ""
        self.suffix: int = 0

    def add(self, record: pysam.AlignedSegment):
        '''add a pysam alignment record to the read cloud, keeping only the relevant info'''
        if self.valid:
            self.positions.append([record.reference_start, record.reference_end])
            self.inserts.append(insert_size(record))
            self.barcode = record.get_tag("BX")
        self.bp.append(record.query_alignment_length)
        self.count.append(record.is_read1 or not record.is_paired)
        self.chromosome = record.reference_name

    def deconvolve(self, cutoff):
        '''
        Process the read cloud and deconvolute using `cutoff` (if it's >0). Deconvolution
        appends `-N` to the barcodes where `N` is an integer (e.g. `-1`, `-2`). Writes to
        stdout and resets the ReadCloud to retain only `barcode`, `suffix`, and `valid`.
        Returns without writing if the Readcloud is empty.
        '''
        if not self.positions:
            return
        result = ""
        # sort alignment extrema by leftmost position and sort the subsequent info the same way
        sort_values = [sublist[0] for sublist in self.positions]
        sorted_indices = np.argsort(sort_values)
        self.positions = [self.positions[i] for i in sorted_indices]
        self.bp = [self.bp[i] for i in sorted_indices]
        self.inserts = [self.inserts[i] for i in sorted_indices]
        self.count = [self.count[i] for i in sorted_indices]

        # instantiate with the first value
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
                BC = f"{self.barcode}-{self.suffix}" if self.suffix > 0 else self.barcode
                result += f"{self.chromosome}\t{BC}\t" + self.stats(start, end, insert, bp, count)

                # start new molecule with current alignment
                self.suffix += 1
                start = self.positions[idx][0]
                end = curr_end
                insert = self.inserts[idx]
                bp = self.bp[idx]
                count = int(self.count[idx])

        # write final molecule
        BC = f"{self.barcode}-{self.suffix}" if self.suffix > 0 else self.barcode
        result += f"{self.chromosome}\t{BC}\t" + self.stats(start, end, insert, bp, count)
        sys.stdout.write(result)
        self.reset()

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
            cov_ins = max(0, round(min(insert / inferred, 1.0), 5))
        except ZeroDivisionError:
            cov_bp = 0
            cov_ins = 0
        return f"{max(1,count)}\t{start}\t{end}\t{inferred}\t{bp}\t{insert}\t{cov_bp}\t{cov_ins}\n"

    def reset(self):
        '''Reset the values in the class, keeping only `valid`, `barcode` and `suffix`'''
        self.chromosome = ""
        self.positions = []
        self.bp = []
        self.inserts = []
        self.count = []

def writestats(x: dict[str,ReadCloud], thresh):
    '''write to file the bx stats dictionary as a table'''
    for cloud in x.values():
        if not cloud.valid:
            sys.stdout.write(f"{cloud.chromosome}\tinvalid\t" + cloud.stats(0, 0, 0, sum(cloud.bp), sum(cloud.count)))
            cloud.reset()
            continue
        cloud.deconvolve(thresh)

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

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.option('-d', '--distance-threshold', default = 0, show_default = True, type = click.IntRange(min = 0, max_open=True), help = 'Calculate statistics assuming this distance threshold for linking alignments sharing a barcode')
@click.argument('input', required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def bx_stats_sam(distance_threshold, input):
    """
    Linked-read metrics from alignment files

    The alignment file is expected to be in "standard" linked-read format, that is,
    the barcode is contained in the BX:Z tag and the barcode validation is stored
    in the VX:i tag as 0 (invalid) or 1 (valid). Metrics include (per molecule): 
    number of reads, position start, position end, length of molecule inferred from
    alignments, total aligned basepairs, total, length of inferred inserts, molecule
    coverage (%) based on aligned bases, molecule coverage (%) based on total inferred
    insert length. Input file *must be coordinate sorted*.
    """
    sys.stdout.write("contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n")
    with pysam.AlignmentFile(input, require_index=False) as alnfile:
        d = {}
        LAST_CONTIG = None

        for read in alnfile.fetch(until_eof=True):
            chrom = read.reference_name
            # check if the current chromosome is different from the previous one
            # if so, process the dict
            if LAST_CONTIG and chrom != LAST_CONTIG:
                writestats(d, distance_threshold)
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
        writestats(d, distance_threshold)
