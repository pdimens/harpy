"""Harpy sequence adapter trimming and quality control"""

import os
import rich_click as click
from harpy.validation.fastq import FASTQ
from harpy.common.cli_filetypes import HPCProfile, FASTQfile
from harpy.common.cli_params import FastpParams, SnakemakeParams
from harpy.common.file_ops import filepath
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
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup', panel = "Workflow Options",  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.help_option('--help', hidden = True)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def qc(inputs, output_dir, unlinked, min_length, max_length, trim_adapters, deduplicate, extra_params, threads, snakemake, skip_reports, quiet, hpc, clean, container, setup):
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
    workflow = Workflow("qc", "qc.smk", output_dir, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.notebook_files = ["qc_bx_stats.ipynb"]
    workflow.conda = ["qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc = not unlinked, quiet = quiet > 0)
    if trim_adapters:
        if trim_adapters != "auto":
            if not os.path.isfile(trim_adapters):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not exist. Please check the spelling or verify the location of the file.")
            if not os.access(trim_adapters, os.R_OK):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not have read permissions. Please modify the persmissions of the file to grant read access.")
            trim_adapters = filepath(trim_adapters)
    else:
        trim_adapters = False

    workflow.notebooks["skip"] = skip_reports
    workflow.linkedreads["type"] = fastq.lr_type
    workflow.input(fastq.files)
    workflow.param(trim_adapters, "trim-adapters")
    workflow.param(deduplicate, "deduplicate")
    workflow.param(min_length, "min-len")
    workflow.param(max_length, "max-len")
    if extra_params:
        workflow.param(extra_params, "extra")

    treatment = ", ".join(i for i,j in zip(["deduplicate", "trim adapters"], [deduplicate, trim_adapters]) if j)

    workflow.info = {
        "Samples": fastq.count,
        "Linked-Read Type": fastq.lr_type,
        **({"Treatment" : treatment} if treatment else {"Treatment": "None"}),
        "Output Folder": os.path.relpath(output_dir) + "/",
    }

    workflow.initialize(setup)
