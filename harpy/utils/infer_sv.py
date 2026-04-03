import sys
import click

conversions = {
    "+-": "deletion",
    "--": "inversion",
    "++": "inversion",
    "-+": "duplication"
}

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.option('-f', '--fail', type=click.File(mode='w', lazy = True))
@click.argument('bedfile', required = True, type=click.File())
@click.help_option('--help', hidden = True)
def infer_sv(bedfile, fail):
    """
    Infer variant types from NAIBR bedpe output

    Removes variants with FAIL flags, use optional -f argument to output FAIL variants to a separate file.
    Use --fail to output a file with variants that failed to pass NAIBR filtering thresholds. Writes to stdout.
    """
    # first line, the header
    line = bedfile.readline().strip().split("\t")
    line_corrected = [i.title().replace(" ", "") for i in line]
    line_corrected.append("SV")
    HEADERLINE = "\t".join(line_corrected) + "\n"
    sys.stdout.write(HEADERLINE)
    if fail:
        fail.write(HEADERLINE)
    # the reamining lines
    while True:
        line = bedfile.readline()
        if not line:
            break
        line = line.strip().split("\t")
        inference = conversions[line[6]]
        line.append(inference)
        NEWROW = "\t".join(line) + "\n"
        if "FAIL" not in line:
            sys.stdout.write(NEWROW)
        elif fail:
            fail.write(NEWROW)

