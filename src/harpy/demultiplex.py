import rich_click as click
from .helperfunctions import fetch_file, generate_conda_deps, validate_demuxschema, check_demux_fastq
import subprocess
import shutil
import sys
import os
import re

#try:
#    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
#except:
#    pass

@click.command(no_args_is_help = True)
@click.option('-f', '--file', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'The forward (or reverse) multiplexed FASTQ file')
@click.option('-b', '--samplesheet', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample<tab>barcode')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def gen1(file, samplesheet, threads, snakemake, quiet, print_only):
    """
    Demultiplex Generation 1 haplotagged FASTQ files

    Use one of the four gzipped FASTQ files provided by the sequencer (I1, I2, R1, R2).fastq.gz for the `--file` argument, 
    Harpy will infer the other three. Note: the `--samplesheet` must be tab (or space) delimited and have no header (i.e. no column names).
    """
    check_demux_fastq(file)

    snakefile = fetch_file("demultiplex.gen1.smk")
    inprefix = re.sub(r"[\_\.][IR][12]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", "", os.path.basename(file))
    os.makedirs(f"Demultiplex/{inprefix}/logs/", exist_ok = True)
    # copy2 to keep metadata during copy
    shutil.copy2(snakefile, f"Demultiplex/{inprefix}/logs/demultiplex.gen1.smk")

    validate_demuxschema(samplesheet)

    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile Demultiplex/{inprefix}/logs/demultiplex.gen1.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"infile={file}")
    command.append(f"samplefile={samplesheet}")
    if print_only:
        click.echo(" ".join(command))
    else:
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)