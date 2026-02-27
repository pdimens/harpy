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
        description='Plot an interactive histogram of alignment and/or molecule depths. Outputs an html file of {prefix}.{contig}.depth.html per contig.',
        usage = "plot-depth -p prefix -a depth.bed -m depth.molcov contig1 contig2..."
        )
    parser.add_argument("-m", "--molcov", type=str, help="Molecule coverage file, such as that produced by harpy align or the molecule-coverage script")
    parser.add_argument("-a", "--alncov", type=str, help="Alignment coverage file, such as that produced by harpy align or mosdepth")
    parser.add_argument("-p", "--prefix", default = "sample", type=str, help="Output filename prefix")
    parser.add_argument("contigs", nargs = '*', type = str, help = "Name(s) of contigs to plot, space-separated (default = 30 largest)")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    coverage = args.alncov
    molcov = args.molcov
    contigs = args.contigs
    prefix = args.prefix

    if not molcov and not coverage:
        parser.error(f"At least one of `-a` or `-m` need to be provided\n")
    if molcov and not os.path.isfile(molcov):
        parser.error(f"File {molcov} does not exist.\n")
    if coverage and not os.path.isfile(coverage):
        parser.error(f"File {coverage} does not exist.\n")

    # moved here to reduce import lag when using CLI
    from harpy.report.theme import palette
    from harpy.report.components import depthplot
    import pandas as pd

    _contigs = None
    if coverage:
        tb = pd.read_table(coverage, header = None)
        if len(tb) == 0:
            parser.error(f"{coverage} is empty")
        else:
            tb.columns = ["Contig", "Position", "Position End", "Read Depth"]
            tb['Read Depth'] = tb['Read Depth'].round(2)
        _contigs = set(tb['Contig'])
        if contigs:
            missing = set(contigs) - _contigs
            if missing:
                print(f"Requested contigs were not found in {os.path.relpath(coverage)}:\n{'\n'.join(missing)}", file = sys.stderr)
                sys.exit(1)

    if molcov:
        tbmol = pd.read_table(molcov, sep = '\t', header = None)
        if len(tbmol) == 0:
            parser.error(f"{molcov} is empty")
        else:
            tbmol.columns = ["Contig", "Position", "Position End", "Molecule Depth"]
            tbmol['Molecule Depth'] = tbmol['Molecule Depth'].round(2)
            if not _contigs:
                _contigs = set(tbmol['Contig'])
                if contigs:
                    missing = set(contigs) - _contigs
                    if missing:
                        print(f"Requested contigs were not found in {os.path.relpath(molcov)}:\n{'\n'.join(missing)}", file = sys.stderr)
                        sys.exit(1)

    if coverage and molcov:
        tb = tb.merge(tbmol, on= ["Contig", "Position", "Position End"])
    elif molcov and not coverage:
        tb = tbmol
    if molcov:
        del tbmol

    if not contigs:
        contigs = tb.groupby('Contig')['Position'].max().nlargest(30).index.tolist()

    tb = tb[tb['Contig'].isin(contigs)].reset_index()
    tb['Genomic Interval (bp)'] = [f"{i}-{j}" for i,j in zip(tb['Position'], tb['Position End'])]

    outdir = os.path.dirname(prefix)
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    grouped = tb.groupby('Contig')
    for name, group in grouped:
        depthplot(group, name).save(f'{prefix}.{name}.depth.html')
