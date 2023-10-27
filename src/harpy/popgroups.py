from .helperfunctions import get_samples_from_fastq, print_error, print_solution

import os
import re
import sys
import glob
import subprocess
import rich_click as click
from rich import print
from rich.panel import Panel

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Input folder Path", help = 'Input folder with fastq or bam files')
@click.option('-o', '--output', type=str, default = "samples.groups", metavar = "Output file name", help = 'Output file name')
def popgroup(directory, output):
    """
    Create a sample grouping file

    With this command you can generate a sample grouping file (for variant calling).
    By default, the output file will be named `samples.groups`.
    **Note** that Harpy cannot reliably infer populations from filenames, therefore all samples
    will be assigned to `pop1`. Please modify this file with appropriate population
    designations.
    """
    try:
        samplenames = getnames(directory, '.bam')
    except:
        full_flist = [i for i in glob.iglob(f"{directory}/*") if not os.path.isdir(i)]
        r = re.compile(".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
        full_fqlist = list(filter(r.match, full_flist))
        fqlist = [os.path.basename(i) for i in full_fqlist]
        bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
        if len(fqlist) == 0:
            print_error(f"No FASTQ or BAM files were detected in [bold]{directory}[/bold]")
            print_solution(
                "Check that FASTQ file endings conform to [green].[/green][[green]F[/green][dim]|[/dim][green]R1[/green]][green].[/green][[green]fastq[/green][dim]|[/dim][green]fq[/green]][green].gz[/green]\nCheck that BAM files end in [green].bam[/green]\nRead the documentation for details: https://pdimens.github.io/harpy/haplotagdata/#naming-conventions"
            )
            exit(1)
        samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

    click.echo(str(len(samplenames)) + f" samples detected in {directory}", file = sys.stderr)
    if os.path.exists(output):
        overwrite = input(f"File {output} already exists, overwrite (no|yes)?  ").lower()
        if (overwrite != "yes") or (overwrite != "y"):
            click.echo("Please suggest a different name for the output file")
            exit(0)
    with open(output, "w") as file:
        for i in samplenames:
            _ = file.write(i + '\tpop1\n') 
    click.echo(f'Created sample population grouping file: {output}', file = sys.stderr, color = True)
    click.echo('\nPlease review it, as all samples have been grouped into a single population\n', file = sys.stderr, color = True)
