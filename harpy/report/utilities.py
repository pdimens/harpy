import gzip
import numpy as np
import pandas as pd

def last_line(filename: str):
    '''Returns the last line of a file. Automatically handles gzip if file ends with case-insensitive `.gz`'''
    if filename.lower().endswith(".gz"):
        with gzip.open(filename, 'rb') as f:
            last_line = None
            for line in f:
                last_line = line
            return last_line.decode("utf-8").strip()
    else:
        with open(filename, 'r') as f:
            last_line = None
            for line in f:
                last_line = line
            return last_line.strip()


def nxx(lengths: list[int]|pd.Series, X:int = 50):
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

def binned_histogramBAK(data: pd.Series, bin_size: int|float, normalize: bool = False):
    '''
    Calculates a binned histogram of counts from the input `data['column']` for bins of size `bin_size`
    with columns ['bin','count']. If `normalize=True`, returns a DataFrame with columns ['bin', 'propportion'].
    '''
    col_max = int(data.max())
    # Creates bins [0-500), [500-1000), etc.
    bins = np.arange(0, col_max + bin_size, bin_size)
    if isinstance(bin_size, int):
        labels = [f"{i}" for i in range(0, col_max, bin_size)]
    else:
        labels = [f"{round(x * bin_size,2)}" for x in range(0, int(1 / bin_size) + 1)]
    binned = pd.cut(data, bins=bins, labels=labels[:len(bins)-1], include_lowest=True)
    colname = 'proportion' if normalize else 'count'

    binned_counts = binned.value_counts(normalize = normalize).sort_index()
    return pd.DataFrame({
                'bin': binned_counts.index,
                colname: binned_counts.values
            })

def binned_histogram(data: pd.Series, bin_size: int|float, normalize: bool = False, max_val = 0, precision = 2):
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

def trunc_digits(x: float,y: int) -> float:
  '''Trucate the input float `x` at decimal digit `y` without rounding'''
  return float(f"%.{y}f" % x)