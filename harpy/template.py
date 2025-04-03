"""Harpy module to create a sample grouping file"""

import os
import re
import sys
import glob
from rich import print as rprint
import rich_click as click
from ._printing import print_error, print_solution, print_notice
from .hpc import hpc_generic, hpc_googlebatch, hpc_lsf, hpc_slurm

docstring = {
    "harpy template": [
        {
            "name": "Input Files",
            "commands": ["groupings", "impute"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "HPC Configurations",
            "commands": ["hpc-generic", "hpc-googlebatch", "hpc-lsf", "hpc-slurm"],
            "panel_styles": {"border_style": "green", "subtitle": "run without arguments"}
        },
    ]
}

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def template():
    """
    Create template files and HPC configs for workflows

    Call any one of the sub-commands with `--help` to see its usage.
    - `groupings` creates the optional grouping file for variant calling
    - `impute` create the parameter input for `harpy impute` 
    - `hpc-*` creates snakemake HPC configs to auto-submit cluster jobs
    """

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/snp/#sample-grouping-file")
@click.option('-o', '--output', type=str, default = "samples.groups", help = "Output file name")
@click.argument('inputdir', required=True, type=click.Path(exists=True, file_okay=False))
def groupings(inputdir, output):
    """
    Create a template sample-grouping file

    This command generates a sample grouping file, like the kind optional for variant calling.
    Provide the input fastq/bam directory at the end of the command.
    By default, the output file will be named `samples.groups`.
    **Note** that Harpy cannot reliably infer populations from filenames, therefore all samples
    will be assigned to `pop1`. Please modify this file with appropriate population
    designations.
    """
    try:
        samplenames = set()
        re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)
        for i in os.listdir(inputdir):
            if i.lower().endswith(".bam") or i.lower().endswith(".sam"):
                samplenames.add(re_ext.sub("", os.path.basename(i)))
        if len(samplenames) < 1:
            raise Exception
    except:
       full_flist = [i for i in glob.iglob(f"{inputdir}/*") if not os.path.isdir(i)]
       r = re.compile(r".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
       full_fqlist = list(filter(r.match, full_flist))
       fqlist = [os.path.basename(i) for i in full_fqlist]
       bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
       if len(fqlist) == 0:
           print_error("no files found", f"No [bold]FASTQ[/] or [bold]BAM[/] files were detected in [blue]{inputdir}[/]")
           print_solution(
               "Check that [bold]FASTQ[/] file endings conform to [green].[/][[green]F[/][dim]|[/dim][green]R1[/]][green].[/][[green]fastq[/][dim]|[/dim][green]fq[/]][green].gz[/]" +
               "\nCheck that [bold]BAM[/] files end in [green].bam[/]"+
               "\nRead the documentation for details: https://pdimens.github.io/harpy/haplotagdata/#naming-conventions"
           )
           sys.exit(1)
       samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

    rprint(f"\n[bold]{len(samplenames)}[/] samples detected in [blue]{inputdir}[blue]\n", file = sys.stderr)
    if os.path.exists(output):
        overwrite = input(f"File {output} already exists, overwrite (no|yes)?  ").lower()
        if overwrite not in ["yes", "y"]:
            click.echo("Please suggest a different name for the output file")
            sys.exit(0)
    with open(output, "w", encoding="utf-8") as file:
        for i in samplenames:
            _ = file.write(i + '\tpop1\n')
    print_notice(f"Created sample population grouping file [blue]{output}[/]. Please review it, as all samples have been grouped into a single population")

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/impute/#parameter-file")
@click.option('-o', '--output', type=str, required = True, help = 'Output file name')
def impute(output):
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
        f"Created template imputation parameter file: [blue]{output}[/]\n" +
        "Modify the model parameters as needed, but [yellow bold]do not add/remove columns."
    )

template.add_command(impute)
template.add_command(groupings)
template.add_command(hpc_generic)
template.add_command(hpc_googlebatch)
template.add_command(hpc_lsf)
template.add_command(hpc_slurm)
