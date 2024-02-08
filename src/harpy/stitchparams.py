import os
import sys
import rich_click as click

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/impute/#parameter-file")
@click.option('-o', '--output', type=str, required = True, metavar = "Output file name", help = 'Name of output STITCH parameter file')
def stitchparams(output):
    """
    Create a template STITCH parameter file

    With this command you can create a template STITCH parameter
    file necessary for imputation via `harpy impute`. The resulting
    file will have generic values and should be modified to be appropriate
    for your study system.
    """
    if os.path.exists(output):
        overwrite = input(f"File {output} already exists, overwrite (no|yes)?  ").lower()
        if overwrite not in ["yes", "y"]:
            click.echo("Please suggest a different name for the output file")
            exit(0)
    with open(output, "w") as file:
        _ = file.write('model\tusebx\tbxlimit\tk\ts\tngen\n')
        _ = file.write('diploid\tTRUE\t50000\t10\t5\t50\n')
        _ = file.write('diploid\tTRUE\t50000\t10\t1\t50\n')
        _ = file.write('diploid\tTRUE\t50000\t15\t10\t100')
    click.echo(f"Created STITCH parameter file: {output}", file = sys.stderr)
    click.echo("Modify the model parameters as needed, but " + '\033[1m' + "DO NOT" + '\033[0m' + " add/remove columns", file = sys.stderr, color = True)