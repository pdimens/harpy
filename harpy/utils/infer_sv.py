import sys
import click

conversions = {
    "+-": "deletion",
    "--": "inversion",
    "++": "inversion",
    "-+": "duplication"
}

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.option('-f', '--fail', type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.argument('bedfile', required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def infer_sv(bedfile, failfile):
    """
    Infer variant types from NAIBR bedpe output

    Removes variants with FAIL flags, use optional -f argument to output FAIL variants to a separate file.
    Use --fail to output a file with variants that failed to pass NAIBR filtering thresholds. Writes to stdout.
    """
    if failfile:
        with open(bedfile, "r", encoding="utf-8") as f, open(failfile, "w", encoding="utf-8") as failout:
            # first line, the header
            line = f.readline().strip().split("\t")
            line_corrected = [i.title().replace(" ", "") for i in line]
            line_corrected.append("SV")
            HEADERLINE = "\t".join(line_corrected) + "\n"
            sys.stdout.write(HEADERLINE)
            failout.write(HEADERLINE)
            # the reamining lines
            while True:
                line = f.readline()
                if not line:
                    break
                line = line.strip().split("\t")
                inference = conversions[line[6]]
                line.append(inference)
                NEWROW = "\t".join(line) + "\n"
                if "FAIL" not in line:
                    sys.stdout.write(NEWROW)
                else:
                    failout.write(NEWROW)
    else:
        with open(bedfile, "r", encoding="utf-8") as f:
            # first line, the header
            line = f.readline().strip().split("\t")
            line_corrected = [i.title().replace(" ", "") for i in line]
            line_corrected.append("SV")
            HEADERLINE = "\t".join(line_corrected) + "\n"
            sys.stdout.write(HEADERLINE)
            # the reamining lines
            while True:
                line = f.readline()
                if not line:
                    break
                line = line.strip().split("\t")
                if "FAIL" not in line:
                    inference = conversions[line[6]]
                    line.append(inference)
                    NEWROW = "\t".join(line) + "\n"
                    sys.stdout.write(NEWROW)
