import click
import os
import sys
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='altair')

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.option("-m", "--molcov", type=click.Path(exists = True, dir_okay=False, resolve_path=True), help = "molecule coverage file, such as that produced by harpy align or the molecule-coverage script")
@click.option("-c", "--coverage", type=click.Path(exists = True, dir_okay=False, resolve_path=True), help = "alignment coverage file, such as that produced by harpy align or mosdepth")
@click.option("-p", "--prefix", default = "sample", type=str, help="Output filename prefix")
@click.argument("contigs", nargs = -1, type = str)
@click.help_option('--help', hidden = True)
def plot_depth(contigs, prefix, molcov, coverage):
    """
    Plot histograms of alignment and/or molecule depths
    
    Outputs one html file of {prefix}.{contig}.depth.html per contig.
    - contigs: name(s) of contigs to plot, space-separated (default = 30 largest)
    """
    if not molcov and not coverage:
        sys.stderr.write(f"At least one of `-a` or `-m` need to be provided\n")
        sys.exit(1)
    if molcov and not os.path.isfile(molcov):
        sys.stderr.write(f"File {molcov} does not exist.\n")
        sys.exit(1)
    if coverage and not os.path.isfile(coverage):
        sys.stderr.write(f"File {coverage} does not exist.\n")
        sys.exit(1)

    # moved here to reduce import lag when using CLI
    from harpy.report.theme import palette
    from harpy.report.components import depthplot
    import pandas as pd

    _contigs = None
    if coverage:
        tb = pd.read_table(coverage, header = None)
        if len(tb) == 0:
            sys.stderr.write(f"{coverage} is empty\n")
            sys.exit(1)
        else:
            tb.columns = ["Contig", "Position", "Position End", "Read Depth"]
            tb['Read Depth'] = tb['Read Depth'].round(2)
        _contigs = set(tb['Contig'])
        if contigs:
            missing = set(contigs) - _contigs
            if missing:
                sys.stderr.write(f"Requested contigs were not found in {os.path.relpath(coverage)}:\n{'\n'.join(missing)}\n")
                sys.exit(1)

    if molcov:
        tbmol = pd.read_table(molcov, sep = '\t', header = None)
        if len(tbmol) == 0:
            sys.stderr.write(f"{molcov} is empty\n")
            sys.exit(1)
        else:
            tbmol.columns = ["Contig", "Position", "Position End", "Molecule Depth"]
            tbmol['Molecule Depth'] = tbmol['Molecule Depth'].round(2)
            if not _contigs:
                _contigs = set(tbmol['Contig'])
                if contigs:
                    missing = set(contigs) - _contigs
                    if missing:
                        sys.stderr.write(f"Requested contigs were not found in {os.path.relpath(molcov)}:\n{'\n'.join(missing)}\n")
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
