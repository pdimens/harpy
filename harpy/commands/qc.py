"""Harpy sequence adapter trimming and quality control"""

import os
import rich_click as click
from harpy.validation.fastq import FASTQ
from harpy.common.cli_filetypes import HPCProfile, FASTQfile
from harpy.common.cli_types_generic import MultiInt, SnakemakeParams
from harpy.common.cli_types_params import FastpParams
from harpy.common.misc import container_ok, filepath
from harpy.common.printing import workflow_info
from harpy.common.workflow import Workflow

docstring = {
    "harpy qc": [
        {
            "name": "Parameters",
            "options": ["--deconvolve", "--deconvolve-params", "--deduplicate", "--extra-params", "--lr-type", "--min-length", "--max-length", "--trim-adapters"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/qc")
@click.option('-c', '--deconvolve', type = MultiInt(4), help = "`k` `w` `d` `a` QuickDeconvolution parameters, comma-separated")
@click.option('-d', '--deduplicate', is_flag = True, default = False, help = 'Identify and remove PCR duplicates')
@click.option('-x', '--extra-params', type = FastpParams(), help = 'Additional Fastp parameters, in quotes')
@click.option('-L', '--lr-type', type = click.Choice(['none', 'haplotagging', 'stlfr','tellseq'], case_sensitive=False), default = "haplotagging", show_default=True, help = "Linked read type\n[none, haplotagging, stlfr, tellseq]")
@click.option('-M', '--max-length', default = 150, show_default = True, type=click.IntRange(min = 30), help = 'Maximum length to trim sequences down to')
@click.option('-m', '--min-length', default = 30, show_default = True, type=click.IntRange(min = 5), help = 'Discard reads shorter than this length')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "QC", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-a', '--trim-adapters', type = str, help = 'Detect and trim adapters')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def qc(inputs, output_dir, lr_type, min_length, max_length, trim_adapters, deduplicate, deconvolve, extra_params, threads, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Adapter removal and other FASTQ preprocessing

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    Disable linked-read specific parts of the workflow with `-L none`.
    
    **Standard trimming**
    - a sliding window from front to tail
    - poly-G tail removal

    **Optional quality checks**
    - `-a` remove adapters
      - accepts `auto` for automatic detection or a `FASTA file` of adapters to remove
    - `-d` removes optical PCR duplicates
      - recommended to skip at this step in favor of barcode-assisted deduplication after alignment
    - `-c` resolves barcodes shared between unrelated sequences
      - off by default, activated with [4 integers](https://github.com/RolandFaure/QuickDeconvolution?tab=readme-ov-file#usage), separated by commas. `21,40,3,0` would be the QuickDeconvolution defaults
      - use `harpy deconvolve` to perform this task separately
    """
    workflow = Workflow("qc", "qc.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["qc_bx_stats.qmd"]
    workflow.conda = ["qc", "r"]

    ## checks and validations ##
    fastq = FASTQ(inputs)
    if trim_adapters:
        if trim_adapters != "auto":
            if not os.path.exists(trim_adapters):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not exist. Please check the spelling or verify the location of the file.")
            if not os.access(trim_adapters, os.R_OK):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not have read permissions. Please modify the persmissions of the file to grant read access.")
            trim_adapters = filepath(trim_adapters)
    else:
        trim_adapters = False

    workflow.config = {
        "workflow" : workflow.name,
        "linkedread_type" : lr_type.lower(),
        "trim_adapters" : trim_adapters,
        "deduplicate" : deduplicate,
        "min_len" : min_length,
        "max_len" : max_length,
        **({'extra': extra_params} if extra_params else {}),
        **({'deconvolve': {
            "kmer_length" : deconvolve[0],
            "window_size" : deconvolve[1],
            "density" : deconvolve[2],
            "dropout" : deconvolve[3]
        }} if deconvolve else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative
        },
        "conda_environments" : workflow.conda,
        "reports" : {"skip": skip_reports},
        "inputs" : fastq.files
    }

    workflow.start_text = workflow_info(
        ("Samples:", fastq.count),
        ("Trim Adapters:", "yes" if trim_adapters else "no"),
        ("Deduplicate:", "yes" if deduplicate else "no"),
        ("Deconvolve:", "yes" if deconvolve else "no"),
        ("Output Folder:", f"{output_dir}/"),
    )

    workflow.initialize(setup_only)
