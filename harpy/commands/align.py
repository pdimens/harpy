"""Harpy align workflows"""

import os
import rich_click as click
from harpy.common.cli_types_generic import ContigList, InputFile, HPCProfile, SnakemakeParams
from harpy.common.cli_types_params import BwaParams, StrobeAlignParams
from harpy.common.misc import container_ok
from harpy.common.parsers import parse_fastq_inputs
from harpy.common.printing import workflow_info
from harpy.common.validations import check_fasta, fasta_contig_match, fastq_has_bx
from harpy.common.workflow import Workflow

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def align():
    """
    Align sequences to a reference genome

    Provide an additional subcommand `bwa` or `strobe` to get more information on using
    those aligners. Both have comparable performance, but `strobe` is typically faster.
    """

module_docstring = {
    "harpy align": [
        {
            "name": "Commands",
            "commands": ["bwa", "strobe"],
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
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
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
    workflow = Workflow("align_bwa", "align_bwa.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["align_stats.qmd", "align_bxstats.qmd"]
    workflow.conda = ["align", "r", "qc"]

    ## checks and validations ##
    fqlist, sample_count = parse_fastq_inputs(inputs, "INPUTS")
    check_fasta(reference)
    if contigs:
        fasta_contig_match(contigs, reference)
    if ignore_bx:
        is_standardized = False
    else:
        is_standardized = fastq_has_bx(fqlist, threads, quiet)

    workflow.config = {
        "workflow" : workflow.name,
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
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "reference": reference,
            "fastq": fqlist
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:",sample_count),
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

@click.command(no_args_is_help = True, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/strobe/")
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', type = StrobeAlignParams(), help = 'Additional aligner parameters, in quotes')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to pass filtering')
@click.option('-d', '--molecule-distance', default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Align/strobealign", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
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
    workflow = Workflow("align_strobe", "align_strobe.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["align_stats.qmd", "align_bxstats.qmd"]
    workflow.conda = ["align", "r", "qc"]

    ## checks and validations ##
    fqlist, sample_count = parse_fastq_inputs(inputs, "INPUTS")
    check_fasta(reference)
    if contigs:
        fasta_contig_match(contigs, reference)
    if ignore_bx:
        is_standardized = False
    else:
        is_standardized = fastq_has_bx(fqlist, threads, quiet)

    workflow.config = {
        "workflow" : workflow.name,
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
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "reference": reference,
            "fastq": fqlist
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:",sample_count),
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

align.add_command(bwa)
align.add_command(strobe)
