import math
import sys

import click
import numpy as np

from harpy.common.file_ops import safe_read

def n_windows(contig_len: int, windowsize: int) -> int:
    return math.ceil(contig_len / windowsize)

def quantify_overlaps(start: int, end: int, windowsize: int, n_bins: int):
    """
    Yield (bin_idx, bp) pairs for every window the molecule [start, end) touches.
    """
    first = start // windowsize
    last = min(end // windowsize, n_bins - 1)
    for b in range(first, last + 1):
        bin_start = b * windowsize
        bp = min(bin_start + windowsize, end) - max(bin_start, start)
        if bp > 0:
            yield b, bp

def print_depth_counts(contig: str, counts: np.ndarray, windowsize: int, contig_len: int):
    n = len(counts)
    for idx in range(n):
        bin_start = idx * windowsize
        bin_end = min(bin_start + windowsize, contig_len)
        span = bin_end - bin_start
        if span > 0:
            sys.stdout.write(f"{contig}\t{bin_start}\t{bin_end}\t{round(counts[idx] / span, 4)}\n")


@click.command(no_args_is_help=True)
@click.option('-w', '--window', default=10000, show_default=True, type=click.IntRange(min=100, max_open=True), help="Window size (in bp) to sum depths over")
@click.argument('fai', required=True, type=click.File())
@click.argument('statsfile', required=True, type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden=True)
def molecule_coverage(fai, statsfile, window):
    '''Calculate molecule coverage from a barcode stats file'''
    
    contigs = {s[0]: int(s[1]) for line in fai for s in (line.split(),)}

    with safe_read(statsfile) as sf:
        header = sf.readline().rstrip().split()
        try:
            IDX_CONTIG = header.index("contig")
            IDX_START  = header.index("start")
            IDX_END    = header.index("end")
        except ValueError as e:
            sys.stderr.write(f"Error: Missing required column — {e}\n")
            sys.exit(1)

        last_contig = None
        counts = None
        n_bins = 0

        for line in sf:
            if line.startswith("#"):
                continue
            parts = line.split()
            contig = parts[IDX_CONTIG]

            if contig != last_contig:
                if last_contig is not None:
                    print_depth_counts(last_contig, counts, window, contigs[last_contig])
                n_bins = n_windows(contigs[contig], window)
                counts = np.zeros(n_bins, dtype=np.int64)
                last_contig = contig

            start = int(parts[IDX_START])
            end   = int(parts[IDX_END])
            for idx, bp in quantify_overlaps(start, end, window, n_bins):
                counts[idx] += bp

        if last_contig is not None:
            print_depth_counts(last_contig, counts, window, contigs[last_contig])