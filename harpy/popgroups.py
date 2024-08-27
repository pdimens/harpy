"""Harpy module to create a sample grouping file"""

import os
import re
import sys
import glob
import rich_click as click
from ._printing import print_error, print_solution, print_notice


@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/snp/#sample-grouping-file")
@click.option('-o', '--output', type=str, default = "samples.groups", metavar = "Output file name", help = 'Output file name, will overwrite existing')
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
        re_ext = re.compile("\.(bam|sam)$", re.IGNORECASE)
        for i in os.listdir(inputdir):
            if i.lower().endswith(".bam") or i.lower().endswith(".sam"):
                samplenames.add(re_ext.sub("", os.path.basename(i)))
        if len(samplenames) < 1:
            raise Exception
    except:
       full_flist = [i for i in glob.iglob(f"{inputdir}/*") if not os.path.isdir(i)]
       r = re.compile(".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
       full_fqlist = list(filter(r.match, full_flist))
       fqlist = [os.path.basename(i) for i in full_fqlist]
       bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
       if len(fqlist) == 0:
           print_error("no files found", f"No FASTQ or BAM files were detected in [bold]{inputdir}[/bold]")
           print_solution(
               "Check that FASTQ file endings conform to [green].[/green][[green]F[/green][dim]|[/dim][green]R1[/green]][green].[/green][[green]fastq[/green][dim]|[/dim][green]fq[/green]][green].gz[/green]" +
               "\nCheck that BAM files end in [green].bam[/green]"+
               "\nRead the documentation for details: https://pdimens.github.io/harpy/haplotagdata/#naming-conventions"
           )
           sys.exit(1)
       samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

    click.echo(f"{len(samplenames)} samples detected in {inputdir}", file = sys.stderr)
    if os.path.exists(output):
        write_text = f"The file [green]{output}[/green] was overwritten."
    else:
        write_text = f"Created sample population grouping file [green]{output}[/green]."
    with open(output, "w", encoding="utf-8") as file:
        for i in samplenames:
            _ = file.write(i + '\tpop1\n')
    print_notice(write_text + " Please review it, as all samples have been grouped into a single population")
