import heapq
import sys
from dataclasses import dataclass
from typing import Optional

import click
from pysam import AlignedSegment, AlignmentFile


# ── helpers ──────────────────────────────────────────────────────────────────

def insert_size(rec: AlignedSegment) -> int:
    if rec.is_paired:
        return rec.query_alignment_length if rec.is_supplementary else max(0, rec.template_length)
    return max(abs(rec.template_length), rec.infer_query_length())


def format_molecule(chrom: str, barcode: str, start: int, end: int, insert: int, bp: int, count: int) -> str:
    inferred = end - start
    if inferred:
        cov_bp  = max(0.0, round(min(bp     / inferred, 1.0), 5))
        cov_ins = max(0.0, round(min(insert / inferred, 1.0), 5))
    else:
        cov_bp = cov_ins = 0.0
    return (
        f"{chrom}\t{barcode}\t{max(1, count)}\t{start}\t{end}\t"
        f"{inferred}\t{bp}\t{insert}\t{cov_bp}\t{cov_ins}\n"
    )


# ── per-barcode state ─────────────────────────────────────────────────────────

@dataclass(slots=True)
class _Mol:
    """Running totals for the molecule currently being assembled."""
    start:    int
    end:      int       # max reference_end seen (molecule extent)
    last_end: int       # reference_end of most-recently-added read (gap detection)
    bp:       int
    insert:   int
    count:    int


class ReadCloud:
    """
    O(1)-memory barcode state. Exploits coordinate-sort order to detect splits
    on-the-fly, replacing the original list-accumulate-then-sort pattern.
    """
    __slots__ = ("barcode", "chromosome", "suffix", "valid", "cutoff", "_mol",
                 "_inv_bp", "_inv_count")

    def __init__(self, barcode: str, chromosome: str, valid: bool = True, initial_suffix: int = 0, cutoff: int = 0):
        self.barcode    = barcode
        self.chromosome = chromosome
        self.suffix     = initial_suffix
        self.valid      = valid
        self._mol: Optional[_Mol] = None
        self._inv_bp    = 0
        self._inv_count = 0
        self.cutoff     = cutoff

    def add(self, record: AlignedSegment) -> str:
        """
        Incorporate one read. Returns a completed-molecule line on a split,
        otherwise returns "". The gap is measured against the previous read's
        end (not the molecule max-end), matching the original deconvolve() semantics.
        """
        self.chromosome = record.reference_name
        bp  = record.query_alignment_length
        ins = insert_size(record)
        cnt = int(not record.is_paired or record.is_read1)

        if not self.valid:
            self._inv_bp    += bp
            self._inv_count += cnt
            return ""

        start = record.reference_start
        end   = record.reference_end

        if self._mol is None:
            self._mol = _Mol(start, end, end, bp, ins, cnt)
            return ""

        if self.cutoff > 0 and (start - self._mol.last_end) > self.cutoff:
            completed = self._emit()
            self.suffix += 1
            self._mol = _Mol(start, end, end, bp, ins, cnt)
            return completed

        self._mol.end      = max(self._mol.end, end)
        self._mol.last_end = end
        self._mol.bp      += bp
        self._mol.insert  += ins
        self._mol.count   += cnt
        return ""

    def flush(self) -> str:
        """Emit whatever is in progress and reset the molecule slot."""
        if not self.valid:
            if self.chromosome:
                return format_molecule(self.chromosome, "invalid", 0, 0, 0, self._inv_bp, self._inv_count)
            return ""
        if self._mol is None:
            return ""
        line = self._emit()
        self._mol = None
        return line

    @property
    def last_end(self) -> int:
        """Most recent read-end; used by the eviction heap."""
        return self._mol.last_end if self._mol else 0

    def _emit(self) -> str:
        bc = f"{self.barcode}-{self.suffix}" if self.suffix else self.barcode
        m = self._mol
        return format_molecule(self.chromosome, bc, m.start, m.end, m.insert, m.bp, m.count)


# ── command ───────────────────────────────────────────────────────────────────

@click.command(no_args_is_help=True, context_settings={"allow_interspersed_args": False})
@click.option('-d', '--distance-threshold', default=0, show_default=True, type=click.IntRange(min=0, max_open=True), help='Distance threshold for splitting molecules sharing a barcode')
@click.argument('input', required=True, type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden=True)
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
    write = sys.stdout.write   # local binding avoids repeated attribute lookup
    write(
        "contig\tmolecule\treads\tstart\tend\tlength_inferred\t"
        "aligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n"
    )

    with AlignmentFile(input, require_index=False) as alnfile:
        clouds: dict[str, ReadCloud] = {}
        # Persists the suffix counter across early evictions so that if a barcode
        # reappears later on the same contig its molecule numbering continues.
        suffix_next: dict[str, int] = {}
        # Min-heap of (last_end, bx). Entries go stale when a cloud is updated;
        # lazy deletion handles this cheaply.
        evict_heap: list[tuple[int, str]] = []
        last_contig: Optional[str] = None
        last_pos = 0

        def flush_all() -> None:
            for _, cloud in clouds.items():
                line = cloud.flush()
                if line:
                    write(line)
            clouds.clear()
            evict_heap.clear()
            suffix_next.clear()   # per-contig; reset on contig boundary

        for read in alnfile.fetch(until_eof=True):
            if (read.is_unmapped
                    or read.is_duplicate
                    or read.is_secondary
                    or (read.is_supplementary and read.reference_name != read.next_reference_name)
                    or not read.get_blocks()):
                continue

            chrom = read.reference_name
            if last_contig and chrom != last_contig:
                flush_all()
                last_pos = 0
            last_contig = chrom
            pos = read.reference_start
            if pos < last_pos:
                sys.stderr.write(
                    "Error: Input alignment file does not appear to be coordinate-sorted, which will result in incorrect output statistics. "
                    "Please sort the alignments by coordinates (e.g. samtools sort) and try again.\n"
                )
                sys.exit(1)
            last_pos = pos
            # ── early eviction ────────────────────────────────────────────
            # Any barcode whose most-recent read ends more than `cutoff` bp
            # behind the current position cannot gain another read to its
            # current molecule. Safe to emit and reclaim now rather than at
            # contig end.
            if distance_threshold > 0:
                while evict_heap and evict_heap[0][0] + distance_threshold < pos:
                    evicted_end, bx = heapq.heappop(evict_heap)
                    cloud = clouds.get(bx)
                    if cloud is None or cloud.last_end > evicted_end:
                        continue   # stale heap entry — cloud was updated after push
                    line = cloud.flush()
                    if line:
                        write(line)
                    suffix_next[bx] = cloud.suffix + 1
                    del clouds[bx]

            # ── barcode lookup ────────────────────────────────────────────
            try:
                bx = read.get_tag("BX")
                if read.get_tag("VX") == 0:
                    raise KeyError
            except KeyError:
                if "invalid" not in clouds:
                    clouds["invalid"] = ReadCloud("invalid", chrom, valid=False)
                clouds["invalid"].add(read)
                continue

            if bx not in clouds:
                clouds[bx] = ReadCloud(bx, chrom, initial_suffix=suffix_next.get(bx, 0), cutoff = distance_threshold)

            line = clouds[bx].add(read)
            if line:
                write(line)

            if distance_threshold > 0:
                heapq.heappush(evict_heap, (clouds[bx].last_end, bx))

        flush_all()