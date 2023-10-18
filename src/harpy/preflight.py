import rich_click as click
import subprocess
import re
import os
import sys
import glob

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with FASTQ files')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def fastq(directory, threads, snakemake, quiet):
    """
    Run validity checks on haplotagged FASTQ files.

    For FASTQ sequence files, it will check that reads have `BX:Z:` tags, that haplotag
    barcodes are propery formatted (`AxxCxxBxxDxx`) and that the comments in the
    read headers conform to the SAM specification of `TAG:TYPE:VALUE`. This **will not**
    fix your data, but it will report the number of reads that feature errors to help
    you diagnose if file formatting will cause downstream issues.
    """
    flist = [os.path.basename(i) for i in glob.iglob(f"{directory}/*") if not os.path.isdir(i)]
    r = re.compile(".*\.f(?:ast)?q\.gz$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, flist))
    if len(fqlist) == 0:
        click.echo(f"\033[1;33mERROR:\033[00m No fastq files with acceptable names found in {directory}", file = sys.stderr, color = True)
        click.echo("Check that the files conform to [.F. | .R1.][.fastq | .fq].gz", file = sys.stderr)
        click.echo("Read the documentation for details: https://pdimens.github.io/harpy/dataformat/#naming-conventions", file = sys.stderr)
        sys.exit(1)

    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/preflight-fastq.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    directory = directory.rstrip("/^")
    command.append(f"seq_directory={directory}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)


@click.command(no_args_is_help = True)
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with FASTQ files')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def bam(directory, threads, snakemake, quiet):
    """
    Run validity checks on haplotagged BAM files

    Files must end in `.bam` (lowercase). For BAM alignment files, it will check that alignments have BX:Z: tags, that haplotag
    barcodes are properly formatted (`AxxCxxBxxDxx`) and that the filename matches the `@RG ID` tag.
    This **will not** fix your data, but it will report the number of records that feature errors  to help
    you diagnose if file formatting will cause downstream issues.
    """
    flist = [i for i in glob.iglob(f"{directory}/*") if not os.path.isdir(i) and i.lower().endswith(".bam")]
    if len(flist) == 0:
        click.echo(f"\033[1;33mERROR:\033[00m No bam files with acceptable names found in {directory}", file = sys.stderr, color = True)
        click.echo("Check that the file names end with .bam", file = sys.stderr)
        click.echo("Read the documentation for details: https://pdimens.github.io/harpy/dataformat/#naming-conventions", file = sys.stderr)
        sys.exit(1)

    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/preflight-bam.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    directory = directory.rstrip("/^")
    command.append(f"seq_directory={directory}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)