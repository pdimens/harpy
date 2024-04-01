import rich_click as click
from pathlib import Path
from .helperfunctions import fetch_file, generate_conda_deps
from .fileparsers import get_samples_from_fastq, parse_fastq_inputs
from .printfunctions import print_error, print_solution, print_notice, print_onstart
import subprocess
import sys
import os
from time import sleep 

@click.command(no_args_is_help = True, epilog= "read the docs for more information: https://pdimens.github.io/harpy/modules/align/bwa/")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for read mapping')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance threshold to separate molecules')
@click.option('-f', '--quality-filter', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), metavar = "Integer", help = 'Minimum mapping quality to pass filtering')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional aligner parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-o', '--output-dir', type = str, default = "Align/bwa", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def bwa(input, output_dir, genome, threads, extra_params, quality_filter, molecule_distance, snakemake, skipreports, quiet, print_only):
    """
    Align sequences to genome using BWA MEM
 
    Provide the input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    BWA is a fast, robust, and reliable aligner that does not use barcodes when mapping.
    Instead, Harpy post-processes the alignments using the specified `--molecule-distance`
    to assign alignments to unique molecules. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/align-bwa.smk')
    command.append("--configfile")
    command.append(f"{workflowdir}/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit(0)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    sn = parse_fastq_inputs(input, f"{workflowdir}/input")
    samplenames = get_samples_from_fastq(f"{workflowdir}/input")
    fetch_file("align-bwa.smk", f"{workflowdir}/")
    fetch_file("assignMI.py", f"{workflowdir}/scripts/")
    fetch_file("bxStats.py", f"{workflowdir}/scripts/")
    for i in ["BxStats", "GenCov"]:
        fetch_file(f"{i}.Rmd", f"{workflowdir}/report/")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"genomefile: {genome}\n")
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplenames: {samplenames}\n")
        config.write(f"quality: {quality_filter}\n")
        config.write(f"molecule_distance: {molecule_distance}\n")
        config.write(f"skipreports: {skipreports}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"workflow_call: {call_SM}\n")
   
    generate_conda_deps()
    print_onstart(
        f"Samples: {len(samplenames)}\nOutput Directory: {output_dir}",
        "align bwa"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

#####----------ema--------------------

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/align/ema")
@click.option('-p', '--platform', type = click.Choice(['haplotag', '10x'], case_sensitive=False), default = "haplotag", show_default=True, help = "Linked read bead technology\n[haplotag, 10x]")
@click.option('-w', '--whitelist', type = click.Path(exists=True, dir_okay=False), help = "Barcode whitelist file for tellseq/10x")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for read mapping')
@click.option('-b', '--ema-bins', default = 500, show_default = True, type = click.IntRange(1,1000), metavar = "Integer", help="Number of barcode bins")
@click.option('-f', '--quality-filter', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), metavar = "Integer", help = 'Minimum mapping quality to pass filtering')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional aligner parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-o', '--output-dir', type = str, default = "Align/ema", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def ema(input, output_dir, platform, whitelist, genome, threads, ema_bins, skipreports, extra_params, quality_filter, snakemake, quiet, print_only):
    """
    Align sequences to a genome using EMA

    Provide the input fastq files and/or directories at the end of the
    command as individual files/folders, using shell wildcards
    (e.g. `data/axolotl*.fastq.gz`), or both.

    EMA may improve mapping, but it also marks split reads as secondary
    reads, making it less useful for variant calling with leviathan.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/align-ema.smk')
    command.append("--configfile")
    command.append(f"{workflowdir}/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit(0)

    platform = platform.lower()
    # the tellseq stuff isn't impremented yet, but this is a placeholder for that... wishful thinking
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

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    sn = parse_fastq_inputs(input, f"{workflowdir}/input")
    samplenames = get_samples_from_fastq(f"{workflowdir}/input")
    fetch_file("align-ema.smk", f"{workflowdir}/")
    fetch_file("bxStats.py", f"{workflowdir}/scripts/")
    for i in ["EmaCount", "EmaGencov", "BxStats"]:
        fetch_file(f"{i}.Rmd", f"{workflowdir}/report/")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"genomefile: {genome}\n")
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplenames: {samplenames}\n")
        config.write(f"quality: {quality_filter}\n")
        config.write(f"platform: {platform}\n")
        config.write(f"EMA_bins: {ema_bins}\n")
        config.write(f"skipreports: {skipreports}\n")
        if whitelist:
            config.write(f"whitelist: {whitelist}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        f"Samples: {len(samplenames)}\nPlatform: {platform}\nOutput Directory: {output_dir}/",
        "align ema"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

#####----------minimap2--------------------

@click.command(no_args_is_help = True, epilog= "read the docs for more information: https://pdimens.github.io/harpy/modules/align/bwa/")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for read mapping')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance threshold to separate molecules')
@click.option('-f', '--quality-filter', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), metavar = "Integer", help = 'Minimum mapping quality to pass filtering')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional aligner parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-o', '--output-dir', type = str, default = "Align/minimap", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def minimap(input, output_dir, genome, threads, extra_params, quality_filter, molecule_distance, snakemake, skipreports, quiet, print_only):
    """
    Align sequences to genome using Minimap2
 
    Provide the input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    Minimap2 is an ultra-fast aligner comparable to bwa for sequences >100bp. This aligner does 
    not use barcodes when mapping. Instead, Harpy post-processes the alignments using the
    specified `--molecule-distance` to assign alignments to unique molecules. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/align-minimap.smk')
    command.append("--configfile")
    command.append(f"{workflowdir}/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit(0)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    sn = parse_fastq_inputs(input, f"{workflowdir}/input")
    samplenames = get_samples_from_fastq(f"{workflowdir}/input")
    fetch_file("align-minimap.smk", f"{workflowdir}/")
    fetch_file("assignMI.py", f"{workflowdir}/scripts/")
    fetch_file("bxStats.py", f"{workflowdir}/scripts/")
    for i in ["BxStats", "GenCov"]:
        fetch_file(f"{i}.Rmd", f"{workflowdir}/report/")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"genomefile: {genome}\n")
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplenames: {samplenames}\n")
        config.write(f"quality: {quality_filter}\n")
        config.write(f"molecule_distance: {molecule_distance}\n")
        config.write(f"skipreports: {skipreports}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"workflow_call: {call_SM}\n")
   
    generate_conda_deps()
    print_onstart(
        f"Samples: {len(samplenames)}\nOutput Directory: {output_dir}",
        "align minimap"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)