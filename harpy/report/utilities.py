from io import StringIO

import numpy as np
import pandas as pd
import polars as pl
import os
import json
import shutil
import sys

from harpy.common.file_ops import safe_read
from harpy.common.printing import HarpyPrint

class StopExecution(Exception):
    '''An exception type to prematurely end a notebook without it being considered an error'''
    def _render_traceback_(self):
        return []

def extract_metric(x: list[str], param: str):
    '''Convenience function to find the relevnant sections of the bcftools.stats file and return a table'''
    selectiontext = "".join(s for s in x if s.startswith(f"{param}\t"))
    if not selectiontext:
        return pd.DataFrame()
    try:
        return pd.read_table(StringIO(selectiontext), sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

def last_line(filename: str) -> str:
    '''Returns the last line of a file. Automatically handles gzip if file ends with case-insensitive `.gz`'''
    with safe_read(filename) as f:
        last_line = None
        for line in f:
            last_line = line
        return last_line.strip()


# this isn't a real thing, just an idea
#def mxx(reads, X: int = 50) -> int:
#    '''
#    Calculate and return the MX value of a list of numbers, where `X` is the kind of MX
#    value you want. MX is the number of molecules containing X percent of all your reads.
#    Sort of like an NX, but specific to linked reads to get an idea of data partitioning.
#    '''
#    threshold = sum(reads) * (X / 100)
#    if isinstance(reads, (pd.Series, pl.Series)):
#        _l = reads.to_list()
#    else:
#        _l = list(reads)
#    _l.sort(reverse=True)
#    cum_sum = 0
#    for j,i in enumerate(_l,1):
#        cum_sum += i
#        if cum_sum >= threshold:
#            return j
#    return len(reads)


def nxx(lengths: list[int]|pd.Series, X:int = 50) -> int:
    '''
    Calculte and return the NX value of a list of numbers, where `X` is
    the kind of NX value you want. For example, `X=50` would return the `N50`.
    '''
    threshold = sum(lengths) * (X/100)
    if isinstance(lengths, pd.Series):
        _l = list(lengths)
    else:
        _l = lengths
    _l.sort(reverse = True)
    cum_sum = 0
    for i in _l:
        cum_sum += i
        if cum_sum >= threshold:
            return i
    return max(lengths)

def nxx_polars(lengths: list[int] | pd.Series | pl.Series, X: int = 50) -> int:
    '''
    Calculate and return the NX value of a list of numbers, where `X` is
    the kind of NX value you want. For example, `X=50` would return the `N50`.
    '''
    threshold = sum(lengths) * (X / 100)
    if isinstance(lengths, (pd.Series, pl.Series)):
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

def binned_histogram(data: pd.Series, bin_size: int|float, normalize: bool = False, max_val = 0, precision = 2) -> pd.DataFrame:
    '''
    Calculates a binned histogram of counts from the input `data['column']` for bins of size `bin_size`
    with columns ['bin','interval','count']. If `normalize=True`, returns a DataFrame with columns ['bin','interval', 'proportion'].
    '''
    col_max = max_val if max_val else data.max().astype(int)
    bins = np.arange(0, col_max + (3*bin_size), bin_size).round(precision)
    labels = []
    for i in bins:
        _i = round(i, precision)
        labels.append(f"{_i}-{round(_i + bin_size, precision)}")
    binned = pd.cut(data, bins=bins, labels=bins.astype(str)[:-1], include_lowest=True)
    colname = 'proportion' if normalize else 'count'
    binned_counts = binned.value_counts(normalize = normalize).sort_index()
    return pd.DataFrame({
                'bin': binned_counts.index,
                'interval': labels[:-1],
                colname: binned_counts.values
            })

def binned_histogram_polars(data: pl.Series, bin_size: int|float, normalize: bool = False, max_val: int|float|None = None, precision = 2) -> pl.DataFrame:
    '''
    Calculates a binned histogram of counts from the input `data` for bins of size `bin_size`
    with columns ['bin','interval','count']. If `normalize=True`, returns a DataFrame with columns ['bin','interval', 'proportion'].
    '''
    col_max = int(data.max()) if max_val is None else max_val
    bins = np.arange(0, col_max + bin_size, bin_size).round(precision)
    #bins = np.arange(0, col_max + (3 * bin_size), bin_size).round(precision)

    labels = [f"{round(i, precision)}-{round(i + bin_size, precision)}" for i in bins]

    # Cut into bins using searchsorted
    bin_indices = np.searchsorted(bins, data.to_numpy(), side='left') - 1
    bin_indices = np.clip(bin_indices, 0, len(bins) - 2)

    colname = 'proportion' if normalize else 'count'

    counts = np.bincount(bin_indices, minlength=len(bins) - 1).astype(float)
    values = counts / counts.sum() if normalize else counts

    return pl.DataFrame({
        'bin': bins[:-1].astype(str),
        'interval': labels[:-1],
        colname: values
    })

def process_variants(df, bin_size=50) -> pd.DataFrame:
    """
    Group variants by binning positions into windows
    """
    # Create binned positions
    df['start_bin'] = (df['Start'] // bin_size) * bin_size
    df['end_bin'] = (df['End'] // bin_size) * bin_size

    # Group by contig, type, and binned positions
    grouped = df.groupby(['Contig', 'Type', 'start_bin', 'end_bin']).agg(
        Start=('Start', 'median'),
        End=('End', 'median'),
        n_samples=('Sample', 'count'),
        Samples=('Sample', list)
    ).reset_index(drop=False)

    # Clean up
    grouped['Start'] = grouped['Start'].astype(int)
    grouped['End'] = grouped['End'].astype(int)
    grouped = grouped[['Contig', 'Start', 'End', 'Type', 'n_samples', 'Samples']]
    grouped.columns = ['Contig', 'Start', 'End', 'Type', 'N Samples', 'Samples']

    return grouped

def trunc_digits(x: float,y: int) -> float:
  '''Trucate the input float `x` at decimal digit `y` without rounding'''
  return float(f"%.{y}f" % x)

def human_format(num):
        if num >= 1e9:
            return f"{num/1e9:.2f}G"
        elif num >= 1e6:
            return f"{num/1e6:.2f}M"
        elif num >= 1e3:
            return f"{num/1e3:.2f}K"
        else:
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


def check_tool(name: str, hint: str) -> None:
    if shutil.which(name) is None:
        hp = HarpyPrint()
        hp.error(
            "Missing dependency",
            f"{name} is not found on the PATH and required to proceed.",
            hint
            )
 