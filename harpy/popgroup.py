"""Harpy module to create a sample grouping file"""

import os
import re
import sys
import glob
from rich import print as rprint
import rich_click as click
from ._printing import print_error, print_solution, print_notice


@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/snp/#sample-grouping-file")
@click.option('-o', '--output', type=str, default = "samples.groups", help = "Output file name")
@click.argument('inputdir', required=True, type=click.Path(exists=True, file_okay=False))
def popgroup(inputdir, output):
    """
    Create a template grouping file for samples

    With this command you can generate a sample grouping file (for variant calling). Provide the 
    input fastq/bam directory at the end of the command.
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
           print_error("no files found", f"No [bold]FASTQ[/bold] or [bold]BAM[/bold] files were detected in [blue]{inputdir}[/blue]")
           print_solution(
               "Check that [bold]FASTQ[/bold] file endings conform to [green].[/green][[green]F[/green][dim]|[/dim][green]R1[/green]][green].[/green][[green]fastq[/green][dim]|[/dim][green]fq[/green]][green].gz[/green]" +
               "\nCheck that [bold]BAM[/bold] files end in [green].bam[/green]"+
               "\nRead the documentation for details: https://pdimens.github.io/harpy/haplotagdata/#naming-conventions"
           )
           sys.exit(1)
       samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

    rprint(f"\n[bold]{len(samplenames)}[/bold] samples detected in [blue]{inputdir}[blue]\n", file = sys.stderr)
    if os.path.exists(output):
        overwrite = input(f"File {output} already exists, overwrite (no|yes)?  ").lower()
        if overwrite not in ["yes", "y"]:
            click.echo("Please suggest a different name for the output file")
            sys.exit(0)
    with open(output, "w", encoding="utf-8") as file:
        for i in samplenames:
            _ = file.write(i + '\tpop1\n')
    print_notice(f"Created sample population grouping file [blue]{output}[/blue]. Please review it, as all samples have been grouped into a single population")
