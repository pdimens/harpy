"""Harpy sequence adapter trimming and quality control"""

import os
import rich_click as click
from harpy.validation.fastq import FASTQ
from harpy.common.cli_filetypes import HPCProfile, FASTQfile
from harpy.common.cli_types_generic import SnakemakeParams
from harpy.common.cli_types_params import FastpParams
from harpy.common.file_ops import filepath
from harpy.common.printing import workflow_info
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/qc")
@click.option('-d', '--deduplicate', panel = "Parameters", is_flag = True, default = False, help = 'Identify and remove PCR duplicates')
@click.option('-x', '--extra-params', panel = "Parameters", type = FastpParams(), help = 'Additional Fastp parameters, in quotes')
@click.option('-M', '--max-length', panel = "Parameters", default = 150, show_default = True, type=click.IntRange(min = 30), help = 'Maximum length to trim sequences down to')
@click.option('-m', '--min-length', panel = "Parameters", default = 30, show_default = True, type=click.IntRange(min = 5), help = 'Discard reads shorter than this length')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "QC", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-a', '--trim-adapters', panel = "Parameters", type = str, help = 'Detect and trim adapters')
@click.option('-U','--unlinked', panel = "Parameters", is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only', panel = "Workflow Options",  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def qc(inputs, output_dir, unlinked, min_length, max_length, trim_adapters, deduplicate, extra_params, threads, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    FASTQ adapter removal, quality filtering, etc.

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    Linked-read presence and type is auto-detected, but you may use `-U` to disable the parts
    of the workflow specific to linked-read data.
    
    **Standard trimming**
    - a sliding window from front to tail
    - poly-G tail removal

    **Optional quality checks**
    - `-a` remove adapters
      - accepts `auto` for automatic detection or a `FASTA file` of adapters to remove
    - `-d` removes optical PCR duplicates
      - recommended to skip at this step in favor of barcode-assisted deduplication after alignment
    """
    workflow = Workflow("qc", "qc.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["qc_bx_stats.qmd"]
    workflow.conda = ["qc", "report"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc = not unlinked)
    if trim_adapters:
        if trim_adapters != "auto":
            if not os.path.isfile(trim_adapters):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not exist. Please check the spelling or verify the location of the file.")
            if not os.access(trim_adapters, os.R_OK):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not have read permissions. Please modify the persmissions of the file to grant read access.")
            trim_adapters = filepath(trim_adapters)
    else:
        trim_adapters = False

    workflow.inputs = fastq.files
    workflow.config = {
        "workflow" : workflow.name,
        "linkedreads" : {
            "type" : fastq.lr_type
        },
        "trim_adapters" : trim_adapters,
        "deduplicate" : deduplicate,
        "min_len" : min_length,
        "max_len" : max_length,
        **({'extra': extra_params} if extra_params else {}),
        "reports" : {"skip": skip_reports},
    }

    workflow.start_text = workflow_info(
        ("Samples:", fastq.count),
        ("Linked-Read Type:", fastq.lr_type),
        ("Trim Adapters:", "yes" if trim_adapters else "no"),
        ("Deduplicate:", "yes" if deduplicate else "no"),
        ("Output Folder:", os.path.relpath(output_dir) + "/"),
    )

    workflow.initialize(setup_only)
