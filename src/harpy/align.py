import rich_click as click
from pathlib import Path
from .helperfunctions import fetch_file, generate_conda_deps, get_samples_from_fastq
from .helperfunctions import print_error, print_solution, print_notice
import subprocess
import sys
import os
from time import sleep 

@click.command(no_args_is_help = True)
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for read mapping')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with sample sequences')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance threshold to separate molecules')
@click.option('-f', '--quality-filter', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), metavar = "Integer", help = 'Minimum mapping quality to pass filtering')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional aligner parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def bwa(genome, threads, directory, extra_params, quality_filter, molecule_distance, snakemake, quiet, print_only):
    """
    Align sequences to genome using BWA MEM
 
    BWA is a fast, reliable, and robust aligner that does not use barcodes to improve mapping.
    Instead, Harpy post-processes the alignments using the specified `--molecule-distance`
    to assign alignments to unique molecules.
    """
    fetch_file("align-bwa.smk", "Align/bwa/workflow/")
    for i in ["BxStats", "BwaGencov"]:
        fetch_file(f"report{i}.Rmd", "Align/bwa/workflow/report/")
    samplenames = get_samples_from_fastq(directory)
    directory = directory.rstrip("/^")
    command = f'snakemake --rerun-incomplete --nolock --use-conda --cores {threads} --directory . --snakefile Align/bwa/workflow/align-bwa.smk'.split()
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
    if print_only:
        click.echo(" ".join(command))
    else:
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)

#####----------ema--------------------

@click.command(no_args_is_help = True)
@click.option('-p', '--platform', type = click.Choice(['haplotag', 'tellseq', '10x'], case_sensitive=False), default = "haplotag", show_default=True, help = "Linked read bead technology\n[haplotag, tellseq, 10x]")
@click.option('-w', '--whitelist', type = click.Path(exists=True, dir_okay=False), help = "Barcode whitelist file for tellseq/10x")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for read mapping')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with sample sequences')
@click.option('-b', '--ema-bins', default = 500, show_default = True, type = click.IntRange(1,1000), metavar = "Integer", help="Number of barcode bins")
@click.option('-f', '--quality-filter', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), metavar = "Integer", help = 'Minimum mapping quality to pass filtering')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional aligner parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def ema(platform, whitelist, genome, threads, ema_bins, directory, extra_params, quality_filter, snakemake, quiet, print_only):
    """
    Align sequences to a genome using EMA

    EMA may improve mapping, but it also marks split reads as secondary
    reads, making it less useful for leviathan variant calling.
    """
    fetch_file("align-ema.smk", "Align/ema/workflow/")
    for i in ["EmaCount", "EmaGencov", "BxStats"]:
        fetch_file(f"report{i}.Rmd", "Align/ema/workflow/report/")

    platform = platform.lower()
    if platform in ["tellseq", "10x"] and not whitelist:
        print_error(f"{platform} technology requires the use of a barcode whitelist.")
        if platform == "10x":
            print_solution("Running EMA requires 10X barcodes provided to [green]--whitelist[/green]. A standard 10X barcode whitelist can be downloaded from [dim]https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes[/dim]")
        else:
            print_solution("Running EMA requires TELLseq barcodes provided to [green]--whitelist[/green]. They can be acquired from the TELL-read software [dim]https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/universal-sequencing-tell-seq-data-analysis-pipeline.html[/dim]")
        exit(1)
    if platform == "haplotag" and whitelist:
        print_notice("Haplotag data does not require barcode whitelists and the whitelist provided as input will be ignored.")
        sleep(3)
    samplenames = get_samples_from_fastq(directory)
    directory = directory.rstrip("/^")
    command = f'snakemake --rerun-incomplete --nolock --use-conda --cores {threads} --directory . --snakefile Align/ema/workflow/align-ema.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    if whitelist:
        command.append(f"whitelist={whitelist}")
    command.append(f"platform={platform}")
    command.append(f"genomefile={genome}")
    command.append(f"quality={quality_filter}")
    command.append(f"samplenames={samplenames}")
    command.append(f"EMA_bins={ema_bins}")
    command.append(f"seq_directory={directory}")

    if extra_params is not None:
        command.append(f"extra={extra_params}")
    if print_only:
        click.echo(" ".join(command))
    else:
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)