"""Harpy align workflows"""

import os
import sys
from time import sleep
from pathlib import Path
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_report, fetch_rule, snakemake_log, launch_snakemake
from .fileparsers import parse_fastq_inputs
from .printfunctions import print_error, print_solution, print_notice
from .validations import check_fasta

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def align():
    """
    Align sample sequences to a reference genome

    The available aligners all retain the linked-read barcode information in the
    resulting output, however `ema` is the only aligner to use the barcode information
    to facilitate the aligning process and can be prohibitively slow. The `strobe`
    aligner is the fastest option and is comparable in accuracy (or better) to `bwa` for
    sequences >100bp.

    Provide an additional subcommand `bwa`, `ema`, or `strobe` to get more information on using
    those aligners.
    """

docstring = {
    "harpy align bwa": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--genome", "--keep-unmapped", "--molecule-distance", "--min-quality"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--depth-window", "--hpc", "--output-dir", "--quiet", "--skipreports", "--snakemake", "--threads", "--help"],
        },
    ],
    "harpy align ema": [
        {
            "name": "Parameters",
            "options": ["--ema-bins", "--extra-params", "--genome", "--keep-unmapped", "--platform", "--min-quality", "--whitelist"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--depth-window", "--hpc", "--output-dir", "--quiet", "--skipreports", "--snakemake", "--threads", "--help"],
        },
    ],
    "harpy align strobe": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--genome", "--keep-unmapped", "--molecule-distance", "--min-quality", "--read-length"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--depth-window", "--hpc", "--output-dir", "--quiet", "--skipreports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog= "See the documentation for more information: https://pdimens.github.io/harpy/modules/align/bwa/")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False, readable=True), required = True, help = 'Genome assembly for read mapping')
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = int, help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', type = str, help = 'Additional bwa mem parameters, in quotes')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), help = 'Minimum mapping quality to pass filtering')
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Base-pair distance threshold to separate molecules')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Align/bwa", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--config-only',  is_flag = True, hidden = True, default = False, help = 'Create the config.yaml file and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def bwa(inputs, output_dir, genome, depth_window, threads, keep_unmapped, extra_params, min_quality, molecule_distance, snakemake, skipreports, quiet, hpc, conda, config_only):
    """
    Align sequences to genome using `BWA MEM`
 
    Provide the input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    BWA is a fast, robust, and reliable aligner that does not use barcodes when mapping.
    Harpy will post-processes the alignments using the specified `--molecule-distance`
    to assign alignments to unique molecules. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/align_bwa.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake is not None:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    check_fasta(genome)
    fetch_rule(workflowdir, "align_bwa.smk")
    fetch_report(workflowdir, "align_stats.Rmd")
    fetch_report(workflowdir, "align_bxstats.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "align_bwa")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: align bwa\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"alignment_quality: {min_quality}\n")
        config.write(f"keep_unmapped: {keep_unmapped}\n")
        config.write(f"molecule_distance: {molecule_distance}\n")
        config.write(f"depth_windowsize: {depth_window}\n")
        config.write(f"skip_reports: {skipreports}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome: {Path(genome).resolve()}\n")
        config.write("  fastq:\n")
        for i in fqlist:
            config.write(f"    - {i}\n")

    if config_only:
        sys.exit(0)

    generate_conda_deps()
    start_text =  f"Samples: {sample_count}\nOutput Directory: {output_dir}\nLog: {sm_log}"
    launch_snakemake(command, "align_bwa", start_text, output_dir, sm_log, quiet)

@click.command(no_args_is_help = True, epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/align/ema")
@click.option('-x', '--extra-params', type = str, help = 'Additional ema align parameters, in quotes')
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = int, help = 'Interval size (in bp) for depth stats')
@click.option('-b', '--ema-bins', default = 500, show_default = True, type = click.IntRange(1,1000), help="Number of barcode bins")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False, readable=True), required = True, help = 'Genome assembly for read mapping')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), help = 'Minimum mapping quality to pass filtering')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Align/ema", show_default=True,  help = 'Output directory name')
@click.option('-p', '--platform', type = click.Choice(['haplotag', '10x'], case_sensitive=False), default = "haplotag", show_default=True, help = "Linked read bead technology\n[haplotag, 10x]")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-l', '--whitelist', type = click.Path(exists=True, dir_okay=False), help = "Barcode whitelist file for 10x linked reads")
@click.option('--config-only',  is_flag = True, hidden = True, default = False, help = 'Create the config.yaml file and exit')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def ema(inputs, output_dir, platform, whitelist, genome, depth_window, keep_unmapped, threads, ema_bins, skipreports, extra_params, min_quality, snakemake, quiet, hpc, conda, config_only):
    """
    Align sequences to genome using `EMA`

    Provide the input fastq files and/or directories at the end of the
    command as individual files/folders, using shell wildcards
    (e.g. `data/axolotl*.fastq.gz`), or both.

    EMA may improve mapping, but it also marks split reads as secondary
    reads, making it less useful for variant calling with leviathan. The barcode
    whitelist is a list of barcodes (in nucleotide format, one per line) that lets EMA know what
    sequences at the beginning of the forward reads are known barcodes.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/align_ema.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake is not None:
        command += snakemake

    platform = platform.lower()
    # the tellseq stuff isn't impremented yet, but this is a placeholder for that... wishful thinking
    if platform in ["tellseq", "10x"] and not whitelist:
        print_error(f"{platform} technology requires the use of a barcode whitelist.")
        if platform == "10x":
            print_solution("Running EMA requires 10X barcodes provided to [green]--whitelist[/green]. A standard 10X barcode whitelist can be downloaded from [dim]https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes[/dim]")
        else:
            print_solution("Running EMA requires TELLseq barcodes provided to [green]--whitelist[/green]. They can be acquired from the TELL-read software [dim]https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/universal-sequencing-tell-seq-data-analysis-pipeline.html[/dim]")
        sys.exit(1)
    if platform == "haplotag" and whitelist:
        print_notice("Haplotag data does not require barcode whitelists and the whitelist provided as input will be ignored.")
        sleep(3)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    check_fasta(genome)
    fetch_rule(workflowdir, "align_ema.smk")
    fetch_report(workflowdir, "align_stats.Rmd")
    fetch_report(workflowdir, "align_bxstats.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "align_ema")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: align ema\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"quality: {min_quality}\n")
        config.write(f"keep_unmapped: {keep_unmapped}\n")
        config.write(f"platform: {platform}\n")
        config.write(f"EMA_bins: {ema_bins}\n")
        config.write(f"depth_windowsize: {depth_window}\n")
        config.write(f"skip_reports: {skipreports}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome: {Path(genome).resolve()}\n")
        if whitelist:
            config.write(f"  whitelist: {Path(whitelist).resolve()}\n")
        config.write("  fastq:\n")
        for i in fqlist:
            config.write(f"    - {i}\n")

    if config_only:
        sys.exit(0)

    generate_conda_deps()
    start_text = f"Samples: {sample_count}\nPlatform: {platform}\nOutput Directory: {output_dir}/\nLog: {sm_log}"
    launch_snakemake(command, "align_ema", start_text, output_dir, sm_log, quiet)

@click.command(no_args_is_help = True, epilog= "See the documentation for more information: https://pdimens.github.io/harpy/modules/align/minimap/")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False, readable=True), required = True, help = 'Genome assembly for read mapping')
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = int, help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', type = str, help = 'Additional aligner parameters, in quotes')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(min = 0, max = 40), help = 'Minimum mapping quality to pass filtering')
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Base-pair distance threshold to separate molecules')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Align/strobealign", show_default=True,  help = 'Output directory name')
@click.option('-l', '--read-length', default = "auto", show_default = True, type = click.Choice(["auto", "50", "75", "100", "125", "150", "250", "400"]), help = 'Average read length for creating index')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--config-only',  is_flag = True, hidden = True, default = False, help = 'Create the config.yaml file and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def strobe(inputs, output_dir, genome, read_length, keep_unmapped, depth_window, threads, extra_params, min_quality, molecule_distance, snakemake, skipreports, quiet, hpc, conda, config_only):
    """
    Align sequences to genome using `strobealign`
 
    Provide the input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    strobealign is an ultra-fast aligner comparable to bwa for sequences >100bp and does 
    not use barcodes when mapping, so Harpy will post-processes the alignments using the
    specified `--molecule-distance` to assign alignments to unique molecules. The `--read-length` is
    an *approximate* parameter and should be one of [`auto`, `50`, `75`, `100`, `125`, `150`, `250`, `400`].
    The alignment process will be faster and take up less disk/RAM if you specify an `-r` value that isn't
    `auto`. If your input has adapters removed, then you should expect the read lengths to be <150.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/align_strobealign.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake is not None:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    check_fasta(genome)
    fetch_rule(workflowdir, "align_strobealign.smk")
    fetch_report(workflowdir, "align_stats.Rmd")
    fetch_report(workflowdir, "align_bxstats.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "align_strobe")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: align strobe\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"genomefile: {genome}\n")
        config.write(f"keep_unmapped: {keep_unmapped}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"quality: {min_quality}\n")
        config.write(f"molecule_distance: {molecule_distance}\n")
        config.write(f"average_read_length: {read_length}\n")
        config.write(f"depth_windowsize: {depth_window}\n")
        config.write(f"skip_reports: {skipreports}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome: {Path(genome).resolve()}\n")
        config.write("  fastq:\n")
        for i in fqlist:
            config.write(f"    - {i}\n")

    if config_only:
        sys.exit(0)

    generate_conda_deps()
    start_text =  f"Samples: {sample_count}\nOutput Directory: {output_dir}\nLog: {sm_log}"
    launch_snakemake(command, "align_strobe", start_text, output_dir, sm_log, quiet)

align.add_command(bwa)
align.add_command(ema)
align.add_command(strobe)