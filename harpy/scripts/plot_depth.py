#! /usr/bin/env python
"""Plot interactive histograms of alignment and/or molecule depths"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='altair')

def main():
    parser = argparse.ArgumentParser(
        prog='plot-depth',
        description='Plot an interactive histogram of alignment and/or molecule depths',
        usage = "plot-depth -c <contigs> -a <depth.bed.cov> -m <depth.molcov>"
        )
    parser.add_argument("-m", "--molcov", type=str, help="Molecule coverage file, such as that produced by harpy align or the molecule-coverage script")
    parser.add_argument("-a", "--alncov", type=str, help="Alignment coverage file, such as that produced by harpy align or mosdepth")
    parser.add_argument("contigs", nargs = -1, type = str, help = "Name(s) of contigs to plot, space-separated (default = 30 largest)")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    # moved here to reduce import lag when using CLI
    from harpy.report.theme import palette
    from harpy.report.components import depthplot
    import pandas as pd

    coverage = args.alncov
    molcov = args.molcov

    if not molcov and not coverage:
        parser.error(f"At least one of `-a` or `-m` need to be provided\n")
    if molcov and not os.path.isfile(molcov):
        parser.error(f"File {molcov} does not exist.\n")
    if coverage and not os.path.isfile(coverage):
        parser.error(f"File {coverage} does not exist.\n")

    if coverage:
        tb = pd.read_table(coverage, header = None)
        if len(tb) == 0:
            parser.error(f"{coverage} is empty")
        else:
            tb.columns = ["Contig", "Position", "Position End", "Read Depth"]
            tb['Read Depth'] = tb['Read Depth'].round(2)

    if molcov:
        tbmol = pd.read_table(molcov, sep = '\t', header = None)
        if len(tbmol) == 0:
            parser.error(f"{molcov} is empty")
        else:
            tbmol.columns = ["Contig", "Position", "Position End", "Molecule Depth"]
            tbmol['Molecule Depth'] = tbmol['Molecule Depth'].round(2)

    contigs = args.contigs
    if not args.contigs:
        contigs = tb.groupby('Contig')['Position'].max().nlargest(30).index.tolist()
    else:
        missing = set(contigs) - set(tb['Contig'])
        if missing:
            print(f"Requested contigs were not found in the input files:\n{missing}", file = sys.stderr)
            sys.exit(1)

    if coverage and molcov:
        tb = tb.merge(tbmol, on= ["Contig", "Position", "Position End"])
        del tbmol
    elif molcov and not coverage:
        tb = tbmol
        del tbmol

    tb = tb[tb['Contig'].isin(contigs)].reset_index()
    tb['Genomic Interval (bp)'] = [f"{i}-{j}" for i,j in zip(tb['Position'], tb['Position End'])]

    grouped = tb.groupby('Contig')
    for name, group in grouped:
        depthplot(group, name).save(f'{name}.depth.html')
