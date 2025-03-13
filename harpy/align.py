"""Harpy align workflows"""

import os
import sys
import yaml
import shutil
from pathlib import Path
import rich_click as click
from ._conda import create_conda_recipes
from ._misc import fetch_report, fetch_rule, snakemake_log
from ._cli_types_generic import convert_to_int, ContigList, InputFile, HPCProfile, SnakemakeParams
from ._cli_types_params import BwaParams, EmaParams, StrobeAlignParams
from ._launch import launch_snakemake, SNAKEMAKE_CMD
from ._parsers import parse_fastq_inputs
from ._printing import print_error, print_solution, print_notice, workflow_info
from ._validations import check_fasta, fasta_contig_match, validate_barcodefile

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
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--depth-window", "--hpc", "--ignore-bx", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy align ema": [
        {
            "name": "Parameters",
            "options": [ "--barcode-list", "--ema-bins", "--extra-params", "--fragment-density", "--genome", "--keep-unmapped", "--platform", "--min-quality"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--depth-window", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy align strobe": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--genome", "--keep-unmapped", "--molecule-distance", "--min-quality", "--read-length"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--depth-window", "--hpc", "--ignore-bx", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/bwa/")
@click.option('-g', '--genome', type=InputFile("fasta", gzip_ok = True), required = True, help = 'Genome assembly for read mapping')
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = int, help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', type = BwaParams(), help = 'Additional bwa mem parameters, in quotes')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to pass filtering')
@click.option('-d', '--molecule-distance', default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Align/bwa", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--ignore-bx',  is_flag = True, default = False, help = 'Ignore parts of the workflow specific to linked-read sequences')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def bwa(inputs, output_dir, genome, depth_window, ignore_bx, threads, keep_unmapped, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Align sample sequences to a reference genome using BWA MEM.
    
    This function sets up and executes a Snakemake workflow to perform sequence alignment using BWA MEM.
    It processes input FASTQ files or directories (supporting shell wildcards), validates the reference genome,
    and optionally verifies contig matches. An HPC YAML configuration file can be provided to customize HPC
    submission, and a positive molecule distance triggers post-alignment barcode deconvolution. When setup_only
    is enabled, the workflow configuration is generated and the function exits without launching the workflow.
    
    Args:
        inputs: Input FASTQ files or directories (supports shell wildcards).
        output_dir: Directory for creating workflow outputs and configuration files.
        genome: Path to the reference genome FASTA file.
        depth_window: Window size used for depth calculations.
        ignore_bx: Flag to ignore barcode information during alignment.
        threads: Number of CPU threads to allocate.
        keep_unmapped: Flag indicating whether to retain unmapped reads.
        extra_params: Additional parameters for customizing the workflow.
        min_quality: Minimum mapping quality threshold.
        molecule_distance: Distance threshold for grouping alignments into unique molecules.
        snakemake: Extra command-line arguments appended to the Snakemake command.
        skip_reports: Flag to skip report generation.
        quiet: Flag to run the workflow in quiet mode.
        hpc: Path to an HPC submission YAML configuration file.
        container: Flag indicating whether to use containerized execution.
        contigs: Optional file specifying contig names for validating the reference genome.
        setup_only: If true, prepares the workflow without launching it.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/align_bwa.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    check_fasta(genome)
    if contigs:
        fasta_contig_match(contigs, genome)
    fetch_rule(workflowdir, "align_bwa.smk")
    fetch_report(workflowdir, "align_stats.qmd")
    fetch_report(workflowdir, "align_bxstats.qmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "align_bwa")
    conda_envs = ["align", "r", "qc"]
    configs = {
        "workflow" : "align bwa",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "alignment_quality" : min_quality,
        "keep_unmapped" : keep_unmapped,
        "depth_windowsize" : depth_window,
        "ignore_bx" : ignore_bx,
        "molecule_distance" : molecule_distance,
        **({'extra': extra_params} if extra_params else {}),
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "genome": Path(genome).resolve().as_posix(),
            "fastq": [i.as_posix() for i in fqlist]
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:",sample_count),
        ("Genome:", genome),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "align_bwa", start_text, output_dir, sm_log, quiet, "workflow/align.bwa.summary")

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/align/ema")
@click.option('-x', '--extra-params', type = EmaParams(), help = 'Additional ema align parameters, in quotes')
@click.option('-d', '--fragment-density',  is_flag = True, show_default = True, default = False, help = 'Perform read fragment density optimization')
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = int, help = 'Interval size (in bp) for depth stats')
@click.option('-b', '--ema-bins', default = 500, show_default = True, type = click.IntRange(1,1000, clamp = True), help="Number of barcode bins")
@click.option('-g', '--genome', type=InputFile("fasta", gzip_ok = True), required = True, help = 'Genome assembly for read mapping')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to pass filtering')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Align/ema", show_default=True,  help = 'Output directory name')
@click.option('-p', '--platform', type = click.Choice(['haplotag', '10x'], case_sensitive=False), default = "haplotag", show_default=True, help = "Linked read bead technology\n[haplotag, 10x]")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-l', '--barcode-list', type = click.Path(exists=True, dir_okay=False), help = "File of known barcodes for 10x linked reads")
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def ema(inputs, output_dir, platform, barcode_list, fragment_density, genome, depth_window, keep_unmapped, threads, ema_bins, skip_reports, extra_params, min_quality, snakemake, quiet, hpc, container, contigs, setup_only):
    """
        Align sequences to a reference genome using the EMA algorithm.
    
        Sets up and executes a Snakemake workflow to align input FASTQ data to a
        reference genome using EMA. The function validates input FASTQ files,
        the genome FASTA, and—depending on the sequencing platform—a required barcode
        list (for "10x" and "tellseq" platforms). For "haplotag" data, any provided
        barcode list is ignored. It configures workflow parameters including fragment
        density optimization, depth window, and alignment quality, and it supports
        HPC configuration and containerization. A YAML configuration file is written
        to the workflow directory, and the workflow is launched unless run in setup-only
        mode.
        
        Args:
            inputs: File paths or directories containing FASTQ data; supports shell wildcards.
            output_dir: Directory for workflow artifacts, configuration, and logs.
            platform: Sequencing platform identifier (e.g., "10x", "tellseq", "haplotag").
            barcode_list: Path to a file with known barcodes (one per line); required for some platforms.
            fragment_density: Option for optimizing fragment density.
            genome: Path to the reference genome in FASTA format.
            depth_window: Size of the window for calculating read depth.
            keep_unmapped: Boolean flag to retain unmapped reads.
            threads: Number of computing threads to use.
            ema_bins: EMA-specific binary or flag parameters.
            skip_reports: Flag to skip generating workflow reports.
            extra_params: Additional parameters to customize the Snakemake command.
            min_quality: Minimum quality threshold for alignments.
            snakemake: Extra arguments to append to the Snakemake command.
            quiet: Boolean to suppress non-critical output.
            hpc: Path to the HPC submission YAML configuration file.
            container: Flag indicating whether to use containerization.
            contigs: Optional specification for contig validation against the reference genome.
            setup_only: If True, prepares the workflow configuration without launching it.
        """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/align_ema.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    platform = platform.lower()
    # the tellseq stuff isn't impremented yet, but this is a placeholder for that... wishful thinking
    if platform in ["tellseq", "10x"] and not barcode_list:
        print_error("missing barcode list", f"{platform} technology requires a list of known barcodes.")
        if platform == "10x":
            print_solution("Running EMA requires 10X barcodes provided to [green]--barcode-list[/green]. A standard 10X barcode list can be downloaded from [dim]https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes[/dim]")
        else:
            print_solution("Running EMA requires TELLseq barcodes provided to [green]--barcode-list[/green]. They can be acquired from the TELL-read software [dim]https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/universal-sequencing-tell-seq-data-analysis-pipeline.html[/dim]")
        sys.exit(1)
    if platform == "haplotag" and barcode_list and not quiet:
        print_notice("Haplotag data does not require a barcode list and the file provided to [green]--barcode-list[/green] will be ignored.")

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    check_fasta(genome)
    if contigs:
        fasta_contig_match(contigs, genome)
    if barcode_list:
        validate_barcodefile(barcode_list, False, quiet)
    fetch_rule(workflowdir, "align_ema.smk")
    fetch_report(workflowdir, "align_stats.qmd")
    fetch_report(workflowdir, "align_bxstats.qmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "align_ema")
    conda_envs = ["align", "r", "qc"]

    configs = {
        "workflow" : "align ema",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "alignment_quality" : min_quality,
        "keep_unmapped" : keep_unmapped,
        "fragment_density_optimization": fragment_density,
        "depth_windowsize" : depth_window,
        "platform" : platform,
        "EMA_bins" : ema_bins,
        **({'extra': extra_params} if extra_params else {}),
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "genome": Path(genome).resolve().as_posix(),
            **({'barcode_list': Path(barcode_list).resolve().as_posix()} if barcode_list else {}),
            "fastq": [i.as_posix() for i in fqlist]
        }
    }

    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:",sample_count),
        ("Genome:", genome),
        ("Platform:", platform),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "align_ema", start_text, output_dir, sm_log, quiet, "workflow/align.ema.summary")

@click.command(epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/strobe/")
@click.option('-g', '--genome', type=InputFile("fasta", gzip_ok = True), required = True, help = 'Genome assembly for read mapping')
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = int, help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', type = StrobeAlignParams(), help = 'Additional aligner parameters, in quotes')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to pass filtering')
@click.option('-d', '--molecule-distance', default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Align/strobealign", show_default=True,  help = 'Output directory name')
@click.option('-l', '--read-length', default = "auto", show_default = True, type = click.Choice(["auto", "50", "75", "100", "125", "150", "250", "400"]), help = 'Average read length for creating index')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--ignore-bx',  is_flag = True, default = False, help = 'Ignore parts of the workflow specific to linked-read sequences')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def strobe(inputs, output_dir, genome, read_length, ignore_bx, keep_unmapped, depth_window, threads, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Align sequences using strobealign and configure the Snakemake workflow.
    
    This function sets up and launches a Snakemake workflow for aligning sequencing data with
    strobealign, an ultra-fast aligner optimized for reads longer than 100bp that does not use
    barcodes by default. It processes input FASTQ files or directories (supporting shell wildcards),
    validates the reference genome (and optionally checks contigs), and prepares a YAML configuration
    file that includes parameters such as the approximate read length, alignment quality, and molecule
    distance for barcode deconvolution. If an HPC submission configuration is provided via the `hpc`
    parameter, the corresponding YAML file is copied into a dedicated subdirectory and linked in the
    workflow command. The function also creates necessary conda recipes and, unless run in setup-only
    mode, launches the alignment workflow.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/align_strobealign.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    check_fasta(genome)
    if contigs:
        fasta_contig_match(contigs, genome)
    fetch_rule(workflowdir, "align_strobealign.smk")
    fetch_report(workflowdir, "align_stats.qmd")
    fetch_report(workflowdir, "align_bxstats.qmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "align_strobe")
    conda_envs = ["align", "r", "qc"]
    configs = {
        "workflow" : "align strobe",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "alignment_quality" : min_quality,
        "keep_unmapped" : keep_unmapped,
        "ignore_bx": ignore_bx,
        "average_read_length": read_length,
        "depth_windowsize" : depth_window,
        "molecule_distance" : molecule_distance,
        **({'extra': extra_params} if extra_params else {}),
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "genome": Path(genome).resolve().as_posix(),
            "fastq": [i.as_posix() for i in fqlist]
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:",sample_count),
        ("Genome:", genome),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "align_strobe", start_text, output_dir, sm_log, quiet, "workflow/align.strobealign.summary")

align.add_command(bwa)
align.add_command(ema)
align.add_command(strobe)
