"""Harpy module to create template imputation parameter file"""

import os
import sys
import rich_click as click
from ._printing import print_notice

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/impute/#parameter-file")
@click.option('-o', '--output', type=str, required = True, help = 'Output file name')
def imputeparams(output):
    """
    Create a template imputation parameter file

    With this command you can create a template parameter
    file necessary for imputation via `harpy impute`. The resulting
    file will have generic values and should be modified to be appropriate
    for your study system.
    """
    if os.path.exists(output):
        overwrite = input(f"File {output} already exists, overwrite (no|yes)?  ").lower()
        if overwrite not in ["yes", "y"]:
            click.echo("Please suggest a different name for the output file")
            sys.exit(0)
    with open(output, "w", encoding="utf-8") as file:
        _ = file.write('name\tmodel\tusebx\tbxlimit\tk\ts\tngen\n')
        _ = file.write('k10_ng50\tdiploid\tTRUE\t50000\t10\t1\t50\n')
        _ = file.write('k1_ng30\tdiploid\tTRUE\t50000\t5\t1\t30\n')
        _ = file.write('high_ngen\tdiploid\tTRUE\t50000\t15\t1\t100')
    print_notice(
        f"Created template imputation parameter file: [blue]{output}[/blue]\n" +
        "Modify the model parameters as needed, but [yellow bold]do not add/remove columns."
    )