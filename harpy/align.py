"""Harpy align workflows"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._conda import create_conda_recipes
from ._misc import fetch_report, fetch_rule, instantiate_dir, setup_snakemake, write_workflow_config
from ._cli_types_generic import ContigList, InputFile, HPCProfile, SnakemakeParams
from ._cli_types_params import BwaParams, EmaParams, StrobeAlignParams
from ._launch import launch_snakemake
from ._parsers import parse_fastq_inputs
from ._printing import print_error, print_solution, print_notice, workflow_info
from ._validations import check_fasta, fasta_contig_match, fastq_has_bx, validate_barcodefile

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def align():
    """
    Align sequences to a reference genome

    | aligner     | linked-read aware | speed | best for |
    |:------------|:-----------------:|:-----:|:--------:|
    | bwa         |        🗙         | fair  |  <600bp  |
    | ema         |         ✔         | slow  |  <600bp  |
    | strobealign |        🗙         | fast  |  >100bp  |

    Provide an additional subcommand `bwa`, `ema`, or `strobe` to get more information on using
    those aligners.
    """

module_docstring = {
    "harpy align": [
        {
            "name": "Commands",
            "commands": ["bwa", "ema", "strobe"],
            "panel_styles": {"border_style" : "blue"}
        }
    ]
}

docstring = {
    "harpy align bwa": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--keep-unmapped", "--molecule-distance", "--min-quality"],
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
            "options": [ "--barcode-list", "--ema-bins", "--extra-params", "--fragment-density", "--keep-unmapped", "--platform", "--min-quality"],
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
            "options": ["--extra-params", "--keep-unmapped", "--molecule-distance", "--min-quality", "--read-length"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--depth-window", "--hpc", "--ignore-bx", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/bwa/")
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', type = BwaParams(), help = 'Additional bwa mem parameters, in quotes')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to pass filtering')
@click.option('-d', '--molecule-distance', default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Align/bwa", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--ignore-bx',  is_flag = True, default = False, help = 'Ignore parts of the workflow specific to linked-read sequences')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=InputFile("fasta", gzip_ok = True), required = True, nargs = 1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def bwa(reference, inputs, output_dir, depth_window, ignore_bx, threads, keep_unmapped, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Align sequences to reference genome using BWA MEM
 
    Provide the reference fasta followed by input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    BWA is a fast, robust, and reliable aligner that does not use barcodes when mapping.
    Harpy will post-processes the alignments using the specified `--molecule-distance`
    to assign alignments to unique molecules. Use a value >`0` for `--molecule-distance` to have
    harpy perform alignment-distance based barcode deconvolution.
    """
    workflow = "align_bwa"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    fqlist, sample_count = parse_fastq_inputs(inputs, "INPUTS")
    check_fasta(reference)
    if contigs:
        fasta_contig_match(contigs, reference)
    if ignore_bx:
        is_standardized = False
    else:
        is_standardized = fastq_has_bx(fqlist, threads, quiet)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        workflow,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "align_bwa.smk")
    fetch_report(workflowdir, "align_stats.qmd")
    fetch_report(workflowdir, "align_bxstats.qmd")

    conda_envs = ["align", "r", "qc"]
    configs = {
        "workflow" : workflow,
        "alignment_quality" : min_quality,
        "keep_unmapped" : keep_unmapped,
        "depth_windowsize" : depth_window,
        "barcodes": {
            "ignore" : ignore_bx,
            "standardized": is_standardized,
            "distance_threshold" : molecule_distance,
        },
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "reference": reference,
            "fastq": fqlist
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:",sample_count),
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/align.bwa.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/align/ema")
@click.option('-x', '--extra-params', type = EmaParams(), help = 'Additional ema align parameters, in quotes')
@click.option('-d', '--fragment-density',  is_flag = True, show_default = True, default = False, help = 'Perform read fragment density optimization')
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-b', '--ema-bins', default = 500, show_default = True, type = click.IntRange(1,1000, clamp = True), help="Number of barcode bins")
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to pass filtering')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Align/ema", show_default=True,  help = 'Output directory name')
@click.option('-p', '--platform', type = click.Choice(['haplotagging', '10x'], case_sensitive=False), default = "haplotagging", show_default=True, help = "Linked read type\n[haplotagging, 10x]")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-l', '--barcode-list', type = click.Path(exists=True, dir_okay=False, resolve_path=True), help = "File of known barcodes for 10x linked reads")
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=InputFile("fasta", gzip_ok = True), required = True, nargs = 1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def ema(reference, inputs, output_dir, platform, barcode_list, fragment_density, depth_window, keep_unmapped, threads, ema_bins, skip_reports, extra_params, min_quality, snakemake, quiet, hpc, container, contigs, setup_only):
    """
    Align sequences to reference genome using EMA

    Provide the reference fasta followed by the fastq files and/or directories at the end of the
    command as individual files/folders, using shell wildcards
    (e.g. `data/axolotl*.fastq.gz`), or both.

    EMA may improve mapping, but it also marks split reads as secondary
    reads, making it less useful for variant calling with leviathan. The barcode
    list is a file of known barcodes (in nucleotide format, one per line) that lets EMA know what
    sequences at the beginning of the forward reads are known barcodes.
    """
    workflow = "align_ema"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    platform = platform.lower()
    # the tellseq stuff isn't impremented yet, but this is a placeholder for that (wishful thinking)
    if platform in ["tellseq", "10x"] and not barcode_list:
        print_error("missing barcode list", f"{platform} technology requires a list of known barcodes.")
        if platform == "10x":
            print_solution("Running EMA requires 10X barcodes provided to [green]--barcode-list[/]. A standard 10X barcode list can be downloaded from [dim]https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes[/dim]")
        else:
            print_solution("Running EMA requires TELLseq barcodes provided to [green]--barcode-list[/]. They can be acquired from the TELL-read software [dim]https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/universal-sequencing-tell-seq-data-analysis-pipeline.html[/dim]")
        sys.exit(1)
    if platform == "haplotagging" and barcode_list and not quiet:
        print_notice("Haplotagging data does not require a barcode list and the file provided to [green]--barcode-list[/] will be ignored.")
    fqlist, sample_count = parse_fastq_inputs(inputs, "INPUTS")
    check_fasta(reference)
    if contigs:
        fasta_contig_match(contigs, reference)
    if barcode_list:
        validate_barcodefile(barcode_list, False, quiet, gzip_ok=False, haplotag_only=True)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        workflow,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "align_ema.smk")
    fetch_report(workflowdir, "align_stats.qmd")
    fetch_report(workflowdir, "align_bxstats.qmd")

    conda_envs = ["align", "r", "qc"]
    configs = {
        "workflow" : workflow,
        "alignment_quality" : min_quality,
        "keep_unmapped" : keep_unmapped,
        "fragment_density_optimization": fragment_density,
        "depth_windowsize" : depth_window,
        "platform" : platform,
        "EMA_bins" : ema_bins,
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "reference": reference,
            **({'barcode_list': barcode_list} if barcode_list else {}),
            "fastq": fqlist
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:",sample_count),
        ("Reference:", os.path.basename(reference)),
        ("Platform:", platform),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/align.ema.summary")

@click.command(no_args_is_help = True, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/strobe/")
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', type = StrobeAlignParams(), help = 'Additional aligner parameters, in quotes')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to pass filtering')
@click.option('-d', '--molecule-distance', default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Align/strobealign", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--ignore-bx',  is_flag = True, default = False, help = 'Ignore parts of the workflow specific to linked-read sequences')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=InputFile("fasta", gzip_ok = True), nargs = 1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def strobe(reference, inputs, output_dir, ignore_bx, keep_unmapped, depth_window, threads, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Align sequences to reference genome using strobealign
 
    Provide the reference fasta followed by the input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    strobealign is an ultra-fast aligner comparable to bwa for sequences >100bp and does 
    not use barcodes when mapping. Use a value >`0` for `--molecule-distance` to have
    harpy perform alignment-distance based barcode deconvolution.
    """
    workflow = "align_strobe"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    fqlist, sample_count = parse_fastq_inputs(inputs, "INPUTS")
    check_fasta(reference)
    if contigs:
        fasta_contig_match(contigs, reference)
    if ignore_bx:
        is_standardized = False
    else:
        is_standardized = fastq_has_bx(fqlist, threads, quiet)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        workflow,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "align_strobe.smk")
    fetch_report(workflowdir, "align_stats.qmd")
    fetch_report(workflowdir, "align_bxstats.qmd")

    conda_envs = ["align", "r", "qc"]
    configs = {
        "workflow" : workflow,
        "alignment_quality" : min_quality,
        "keep_unmapped" : keep_unmapped,
        "depth_windowsize" : depth_window,
        "barcodes": {
            "ignore" : ignore_bx,
            "standardized": is_standardized,
            "distance_threshold" : molecule_distance,
        },
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "reference": reference,
            "fastq": fqlist
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:",sample_count),
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/align.strobealign.summary")

align.add_command(bwa)
align.add_command(ema)
align.add_command(strobe)
