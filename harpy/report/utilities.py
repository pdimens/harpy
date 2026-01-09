import gzip
from altair import Data
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
        lengths.sort_values(ascending=False, inplace = True)
    else:
        _ = lengths.sort(reverse = True)
    cum_sum = 0 
    for i in lengths:
        cum_sum += i
        if cum_sum >= threshold:
            return i

def binned_histogram(data: pd.Series, bin_size: int|float, normalize: bool = False):
    '''
    Calculates a binned histogram of counts from the input `data['column']` for bins of size `bin_size`
    with columns ['bin','count']. If `max_val` is provided, returns percents of `max_val` instead of a count
    with columns ['bin', 'propportion'].
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

def trunc_digits(x: float,y: int) -> float:
  '''Trucate the input float `x` at decimal digit `y` without rounding'''
  return float(f"%.{y}f" % x)