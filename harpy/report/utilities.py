from typing import Literal
from io import StringIO

import numpy as np
import polars as pl
import os
import json

from harpy.common.file_ops import safe_read

class StopExecution(Exception):
    '''An exception type to prematurely end a notebook without it being considered an error'''
    def _render_traceback_(self):
        return []

def extract_metric(x: list[str], param: str):
    '''Convenience function to find the relevnant sections of the bcftools.stats file and return a table'''
    selectiontext = "".join(s for s in x if s.startswith(f"{param}\t"))
    if not selectiontext:
        return pl.DataFrame()
    try:
        return pl.read_csv(
            StringIO(selectiontext), separator="\t", has_header=False,
            infer_schema_length=None, null_values=["nan", "-nan"]
        )
    except pl.exceptions.NoDataError:
        return pl.DataFrame()

def last_line(filename: str) -> str:
    '''Returns the last line of a file. Automatically handles gzip if file ends with case-insensitive `.gz`'''
    with safe_read(filename) as f:
        last_line = None
        for line in f:
            last_line = line
        return last_line.strip()

def nxx_polars(lengths: list[int] | pl.Series, X: int = 50) -> int:
    '''
    Calculate and return the NX value of a list of numbers, where `X` is
    the kind of NX value you want. For example, `X=50` would return the `N50`.
    '''
    threshold = sum(lengths) * (X / 100)
    if isinstance(lengths, pl.Series):
        _l = lengths.to_list()
    else:
        _l = list(lengths)
    _l.sort(reverse=True)
    cum_sum = 0
    for i in _l:
        cum_sum += i
        if cum_sum >= threshold:
            return i
    return max(lengths)


def binned_histogram(data: pl.Series, bin_size: int|float, normalize: bool = False, max_val: int|float|None = None, precision = 2) -> pl.DataFrame:
    '''
    Calculates a binned histogram of counts from the input `data` for bins of size `bin_size`
    with columns ['bin','interval','count']. If `normalize=True`, returns a DataFrame with columns ['bin','interval', 'proportion'].
    '''
    col_max: float | int = int(data.max()) if max_val is None else max_val
    bins = np.arange(0, col_max + bin_size, bin_size).round(precision)
    #bins = np.arange(0, col_max + (3 * bin_size), bin_size).round(precision)

    labels: list[str] = [f"{round(i, precision)}-{round(i + bin_size, precision)}" for i in bins]

    # Cut into bins using searchsorted
    bin_indices = np.searchsorted(bins, data.to_numpy(), side='left') - 1
    bin_indices = np.clip(bin_indices, 0, len(bins) - 2)

    colname: Literal['count', 'proportion'] = 'proportion' if normalize else 'count'

    counts = np.bincount(bin_indices, minlength=len(bins) - 1).astype(float)
    values = counts / counts.sum() if normalize else counts

    return pl.DataFrame({
        'bin': bins[:-1].astype(str),
        'interval': labels[:-1],
        colname: values
    })

def process_variants(df: pl.DataFrame, bin_size: int = 50) -> pl.DataFrame:
    """
    Group variants by binning positions into windows
    """
    return (
        df.group_by(
            'Contig', 'Type',
            start_bin=(pl.col('Start') // bin_size) * bin_size,
            end_bin=(pl.col('End') // bin_size) * bin_size,
        )
        .agg(
            pl.col('Start').median().cast(pl.Int64),
            pl.col('End').median().cast(pl.Int64),
            pl.col('Sample').count().alias('N Samples'),
            pl.col('Sample').alias('Samples'),
        )
        .sort('Contig', 'Type', 'start_bin', 'end_bin')
        .select('Contig', 'Start', 'End', 'Type', 'N Samples', 'Samples')
    )

def trunc_digits(x: float,y: int) -> float:
  '''Trucate the input float `x` at decimal digit `y` without rounding'''
  return float(f"%.{y}f" % x)

def human_format(num):
    match num:
        case num if num >= 1e9:
            return f"{num/1e9:.2f}G"
        case num if num >= 1e6:
            return f"{num/1e6:.2f}M"
        case num if num >= 1e3:
            return f"{num/1e3:.2f}K"
        case _:
            return str(num)

class FastpResults():
    def __init__(self, report_dir):
        '''Collect all JSON report files'''
        stats = {
            "Sample"         : [],
            "Reads (Before)" : [],
            "Reads (After)"  : [],
            "Bases (Before)" : [],
            "Bases (After)"  : [],
            "% Q20 (Before)" : [],
            "% Q20 (After)"  : [],
            "% Q30 (Before)" : [],
            "% Q30 (After)"  : [],
            "% GC (Before)"  : [],
            "% GC (After)"   : []
        }
        mean_qual_curves = {"before" : {}, "after" : {}}
        self.max_len = 0
        self.max_len_after = 0
        gc_curves = {"before" : {}, "after" : {}}
        mean_qual_curves_r2 = {"before" : {}, "after" : {}}
        self.max_len_r2 = 0
        self.max_len_r2_after = 0
        gc_curves_r2 = {"before" : {}, "after" : {}}
        json_files = [f for f in os.listdir(report_dir) if f.endswith('.json')]
        self.samples = len(json_files)
        for jf in json_files:
            path = os.path.join(report_dir, jf)
            samplename = jf.replace('.fastp.json', '')
            with open(path) as f:
                data = json.load(f)
                summary = data.get('summary', {})
                before = summary.get('before_filtering', {})
                after = summary.get('after_filtering', {})
                
                # Extract quality and GC curves for read1
                qual_curve_before = data.get('read1_before_filtering', {}).get('quality_curves', {}).get('mean', [])
                qual_curve_after = data.get('read1_after_filtering', {}).get('quality_curves', {}).get('mean', [])
                gc_curve_before = data.get('read1_before_filtering', {}).get('content_curves', {}).get('GC', [])
                gc_curve_after = data.get('read1_after_filtering', {}).get('content_curves', {}).get('GC', [])
                mean_qual_curves["before"][samplename] = qual_curve_before
                mean_qual_curves["after"][samplename] = qual_curve_after
                gc_curves["before"][samplename] = gc_curve_before
                gc_curves["after"][samplename] = gc_curve_after

                self.max_len = max(self.max_len, len(qual_curve_before))
                self.max_len_after = max(len(qual_curve_after), self.max_len_after)

                # Extract quality and GC curves for read2 if present
                qual_curve_before_r2 = data.get('read2_before_filtering', {}).get('quality_curves', {}).get('mean', [])
                qual_curve_after_r2 = data.get('read2_after_filtering', {}).get('quality_curves', {}).get('mean', [])
                gc_curve_before_r2 = data.get('read2_before_filtering', {}).get('content_curves', {}).get('GC', [])
                gc_curve_after_r2 = data.get('read2_after_filtering', {}).get('content_curves', {}).get('GC', [])
                if qual_curve_before_r2 or qual_curve_after_r2 or gc_curve_before_r2 or gc_curve_after_r2:
                    mean_qual_curves_r2["before"][samplename] = qual_curve_before_r2
                    mean_qual_curves_r2["after"][samplename] = qual_curve_after_r2
                    gc_curves_r2["before"][samplename] = gc_curve_before_r2
                    gc_curves_r2["after"][samplename] = gc_curve_after_r2

                    self.max_len_r2 = max(len(qual_curve_before_r2), self.max_len_r2)
                    self.max_len_r2_after = max(len(qual_curve_after_r2), self.max_len_r2_after)

                    stats["Sample"].append(jf.replace('.fastp.json', ''))
                    stats["Reads (Before)"].append(before.get('total_reads', 0))
                    stats["Reads (After)"].append(after.get('total_reads', 0))
                    stats["Bases (Before)"].append(before.get('total_bases', 0))
                    stats["Bases (After)"].append(after.get('total_bases', 0))
                    stats["% Q20 (Before)"].append(round(before.get('q20_rate', 0) * 100, 2))
                    stats["% Q20 (After)"].append(round(after.get('q20_rate', 0) * 100, 2))
                    stats["% Q30 (Before)"].append(round(before.get('q30_rate', 0) * 100, 2))
                    stats["% Q30 (After)"].append(round(after.get('q30_rate', 0) * 100, 2))
                    stats["% GC (Before)"].append(round(before.get('gc_content', 0) * 100, 2))
                    stats["% GC (After)"].append(round(after.get('gc_content', 0) * 100, 2))
        
        self.stats =  pl.DataFrame(stats)
        self.qual_curves = pl.DataFrame(mean_qual_curves["before"])
        self.qual_curves_after = pl.DataFrame(mean_qual_curves["after"])
        self.qual_curves_r2 = pl.DataFrame(mean_qual_curves_r2["before"])
        self.qual_curves_r2_after = pl.DataFrame(mean_qual_curves_r2["after"])
        self.gc_curves = pl.DataFrame(gc_curves["before"])
        self.gc_curves_after = pl.DataFrame(gc_curves["after"])
        self.gc_curves_r2 = pl.DataFrame(gc_curves_r2["before"])
        self.gc_curves_r2_after = pl.DataFrame(gc_curves_r2["after"])
 