import rich_click as click
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
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance threshold to separate molecules')
@click.option('-f', '--quality-filter', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), metavar = "Integer", help = 'Minimum mapping quality to pass filtering')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional aligner parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def bwa(genome, threads, directory, extra_params, quality_filter, molecule_distance, snakemake, quiet):
    """
    Align sequences to genome using BWA MEM
 
    BWA is a fast, reliable, and robust aligner that does not use barcodes to improve mapping.
    Instead, Harpy post-processes the alignments using the specified `--molecule-distance`
    to assign alignments to unique molecules.
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

    samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])
    directory = directory.rstrip("/^")
    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/align-bwa.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"genomefile={genome}")
    command.append(f"quality={quality_filter}")
    command.append(f"samplenames={samplenames}")
    command.append(f"molecule_distance={molecule_distance}")
    command.append(f"seq_directory={directory}")

    if extra_params is not None:
        command.append(f"extra={extra_params}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True)
@click.option('-g', '--genome', type=click.Path(exists=True), required = True, metavar = "File Path", help = 'Genome assembly for read mapping')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with sample sequences')
@click.option('-b', '--ema-bins', default = 500, show_default = True, type = click.IntRange(1,1000), metavar = "Integer", help="Number of barcode bins")
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance threshold to separate molecules')
@click.option('-f', '--quality-filter', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), metavar = "Integer", help = 'Minimum mapping quality to pass filtering')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional aligner parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def ema(genome, threads, ema_bins, directory, extra_params, quality_filter, molecule_distance, snakemake, quiet):
    """
    Align sequences to a genome using EMA

    EMA may improve mapping, but EMA marks split reads as secondary
    reads, making it less useful for leviathan variant calling.
    Note that `--molecule-distance` is for reporting barcode alignment
    information and does not affect mapping.
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

    samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])
    directory = directory.rstrip("/^")
    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/align-ema.smk'.split()
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
    command.append(f"molecule_distance={molecule_distance}")
    command.append(f"seq_directory={directory}")

    if extra_params is not None:
        command.append(f"extra={extra_params}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)