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
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with raw sample sequences')
@click.option('-l', '--max-length', default = 150, show_default = True, type=int, metavar = "Integer", help = 'Maximum length to trim sequences down to')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional Fastp parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def trim(directory, max_length, extra_params, threads, snakemake, quiet):
    """
    Remove adapters and quality trim sequences
    """
    flist = [os.path.basename(i) for i in glob.iglob(f"{directory}/*") if not os.path.isdir(i)]
    r = re.compile(".*\.f(?:ast)?q\.gz$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, flist))
    if len(fqlist) == 0:
        print(f"\033[1;33mERROR:\033[00m No fastq files with acceptable names found in {directory}", file = sys.stderr)
        print("Check that the files conform to [.F. | .R1.][.fastq | .fq].gz", file = sys.stderr)
        print("Read the documentation for details: https://pdimens.github.io/harpy/dataformat/#naming-conventions", file = sys.stderr)
        sys.exit(1)

    command = f'snakemake --rerun-incomplete --cores {threads} --directory . --snakefile {harpypath}/trim.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    directory = directory.rstrip("/^")
    command.append(f"seq_directory={directory}")
    command.append(f"maxlen={max_length}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    subprocess.run(command)