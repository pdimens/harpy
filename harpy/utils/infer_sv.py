import sys

import click

conversions = {
    "+-": "deletion",
    "--": "inversion",
    "++": "inversion",
    "-+": "duplication"
}

@click.command(no_args_is_help = True)
@click.option('-f', '--fail', type=click.File(mode='w', lazy = True))
@click.argument('bedfile', required = True, type=click.File())
@click.help_option('--help', hidden = True)
def infer_sv(bedfile, fail):
    """
    Infer variant types from NAIBR bedpe output

    Removes variants with FAIL flags, use optional -f argument to output FAIL variants to a separate file.
    Use --fail to output a file with variants that failed to pass NAIBR filtering thresholds. Writes to stdout.
    """
    # Read and validate the header before writing anything.  NAIBR's BEDPE
    # output must contain the orientation field used below (index 6 in the
    # current Harpy intermediate schema).
    header = bedfile.readline()
    if not header:
        raise click.ClickException("The NAIBR BEDPE input is empty; expected a header and at least one record.")
    line = header.rstrip("\r\n").split("\t")
    if len(line) < 7:
        raise click.ClickException(
            f"The NAIBR BEDPE header has {len(line)} columns; at least 7 are required to read the orientation field."
        )
    line_corrected = [i.title().replace(" ", "") for i in line]
    line_corrected.append("SV")
    HEADERLINE = "\t".join(line_corrected) + "\n"
    sys.stdout.write(HEADERLINE)
    if fail:
        fail.write(HEADERLINE)
    # Process remaining records while retaining the original columns.
    for line_number, raw_line in enumerate(bedfile, start=2):
        line = raw_line.rstrip("\r\n").split("\t")
        if len(line) < 7:
            raise click.ClickException(
                f"Malformed NAIBR BEDPE record on line {line_number}: found {len(line)} columns; at least 7 are required."
            )
        orientation = line[6]
        if orientation not in conversions:
            raise click.ClickException(
                f"Unsupported NAIBR orientation {orientation!r} on line {line_number}; expected one of {sorted(conversions)}."
            )
        inference = conversions[orientation]
        line.append(inference)
        NEWROW = "\t".join(line) + "\n"
        if "FAIL" not in line:
            sys.stdout.write(NEWROW)
        elif fail:
            fail.write(NEWROW)
