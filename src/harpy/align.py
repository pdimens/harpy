import rich_click as click
#from .harpymisc import sanitize_fastq
from pathlib import Path
import subprocess
import glob
import sys
import os
import re

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-g', '--genome', type=click.Path(exists=True), required = True, metavar = "File Path", help = 'Genome assembly for read mapping')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with sample sequences')
@click.option('-e', '--ema-bins', default = 500, show_default = True, type = click.IntRange(1,1000), metavar = "Integer", help="Number of barcode bins if using EMA")
@click.option('-m', '--method', default = "bwa", show_default = True, type = click.Choice(["bwa", "ema"], case_sensitive = False), metavar = "String", help = "Method for aligning reads")
@click.option('-f', '--quality-filter', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), metavar = "Integer", help = 'Minimum mapping quality to pass filtering')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional aligner parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def align(genome, threads, method, ema_bins, directory, extra_params, quality_filter, snakemake, quiet):
    """
    Align sample sequences to a reference genome

    EMA improves mapping overall, but note that EMA marks split 
    reads as secondary reads, which makes it less useful for leviathan.

    ## methods
    - **bwa**: uses BWA MEM to align reads, retaining BX tags in the alignments
    - **ema**: uses the BX barcode-aware EMA aligner
    """
    full_flist = [i for i in glob.iglob(f"{directory}/*") if not os.path.isdir(i)]
    r = re.compile(".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    full_fqlist = list(filter(r.match, full_flist))
    fqlist = [os.path.basename(i) for i in full_fqlist]
    bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
    if len(fqlist) == 0:
        print(f"\033[1;33mERROR:\033[00m No fastq files with acceptable names found in {directory}", file = sys.stderr)
        print("Check that the files conform to [.F. | .R1.][.fastq | .fq].gz", file = sys.stderr)
        print("Read the documentation for details: https://pdimens.github.io/harpy/dataformat/#naming-conventions", file = sys.stderr)
        sys.exit(1)

    mapper = method
    ## DEPRECATED ##
    # create relative symlinks for input files with standard naming convention
    #linkdir = f"Align/{mapper}/input"
    # find the basenames with this flexible regex
    #samplenames = sanitize_fastq(full_fqlist, linkdir)  
    flist = [os.path.basename(i) for i in glob.iglob(f"{directory}/*") if not os.path.isdir(i)]
    r = re.compile(".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, flist))
    bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
    samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])
    command = f'snakemake --rerun-incomplete --cores {threads} --directory . --snakefile {harpypath}/align-{mapper}.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"genomefile={genome}")
    command.append(f"quality={quality_filter}")
    command.append(f"samplenames={samplenames}")
    command.append(f"EMA_bins={ema_bins}")
    directory = directory.rstrip("/^")
    command.append(f"seq_directory={directory}")

    if extra_params is not None:
        command.append(f"extra={extra_params}")
    subprocess.run(command)