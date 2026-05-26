from io import StringIO
import numpy as np
import pandas as pd
import polars as pl
from harpy.common.file_ops import safe_read

class StopExecution(Exception):
    '''An exception type to prematurely end a notebook without it being considered an error'''
    def _render_traceback_(self):
        """
        Suppress rendering of the traceback in interactive notebook environments.
        
        This IPython/Jupyter hook returns no traceback frames so that the exception appears to terminate execution without displaying a traceback.
        
        Returns:
            list: An empty list indicating no traceback frames should be shown.
        """
        return []

def extract_metric(x: list[str], param: str):
    """
    Extracts contiguous lines that begin with a given parameter prefix from a list of lines and parses them as a tab-separated table.
    
    Parameters:
        x (list[str]): Lines from a bcftools-style stats file or similar tab-delimited text.
        param (str): Parameter prefix to match at line start (matched as f"{param}\\t").
    
    Returns:
        pandas.DataFrame: Parsed table of matching lines; an empty DataFrame if no matching lines are found or parsing yields no data.
    """
    selectiontext = "".join(s for s in x if s.startswith(f"{param}\t"))
    if not selectiontext:
        return pd.DataFrame()
    try:
        return pd.read_table(StringIO(selectiontext), sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

def last_line(filename: str) -> str:
    """
    Return the last line of a text file, handling gzipped files when the filename ends with `.gz` (case-insensitive).
    
    Parameters:
        filename (str): Path to the file to read. If the filename ends with `.gz` (any case), the file will be read as gzip-compressed.
    
    Returns:
        str: The final line of the file with leading and trailing whitespace removed.
    """
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
    """
    Compute the NX statistic (e.g., N50) from a collection of lengths.
    
    Parameters:
        lengths (sequence or pd.Series): Length values to evaluate; a pandas Series is accepted.
        X (int): Percentile for the NX calculation (e.g., 50 for N50).
    
    Returns:
        int: The length value at which the cumulative sum of lengths reaches at least X% of the total; if the threshold is not reached, returns the maximum length.
    """
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

def binned_histogram_polars(data: pl.Series, bin_size: int|float, normalize: bool = False, max_val = 0, precision = 2) -> pl.DataFrame:
    """
    Create a binned histogram from a Polars Series using fixed-size numeric bins.
    
    Parameters:
        data (pl.Series): Numeric values to bin.
        bin_size (int | float): Width of each bin.
        normalize (bool): If True, return proportions instead of raw counts.
        max_val (int | float): Optional upper bound for bins; if zero or falsy, uses the series maximum.
        precision (int): Decimal places to round bin edges and interval labels.
    
    Returns:
        pl.DataFrame: DataFrame with columns:
          - `bin` (str): Left-edge values of bins as strings.
          - `interval` (str): Human-readable interval labels "left-right".
          - `count` or `proportion` (float): Counts per bin or normalized proportions when `normalize=True`.
    """
    col_max = max_val if max_val else int(data.max())
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
