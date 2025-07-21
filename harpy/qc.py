"""Harpy sequence adapter trimming and quality control"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_report, fetch_rule, instantiate_dir, setup_snakemake, write_workflow_config
from ._cli_types_generic import HPCProfile, MultiInt, SnakemakeParams
from ._cli_types_params import FastpParams
from ._misc import filepath
from ._parsers import parse_fastq_inputs
from ._printing import workflow_info
from ._validations import check_fasta

docstring = {
    "harpy qc": [
        {
            "name": "Parameters",
            "options": ["--deconvolve", "--deconvolve-params", "--deduplicate", "--extra-params", "--min-length", "--max-length", "--trim-adapters"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--ignore-bx", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/qc")
@click.option('-c', '--deconvolve', type = MultiInt(4), help = "`k` `w` `d` `a` QuickDeconvolution parameters, comma-separated")
@click.option('-d', '--deduplicate', is_flag = True, default = False, help = 'Identify and remove PCR duplicates')
@click.option('-x', '--extra-params', type = FastpParams(), help = 'Additional Fastp parameters, in quotes')
@click.option('-M', '--max-length', default = 150, show_default = True, type=click.IntRange(min = 30), help = 'Maximum length to trim sequences down to')
@click.option('-m', '--min-length', default = 30, show_default = True, type=click.IntRange(min = 5), help = 'Discard reads shorter than this length')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "QC", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-a', '--trim-adapters', type = str, help = 'Detect and trim adapters')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--ignore-bx',  is_flag = True, default = False, help = 'Ignore parts of the workflow specific to linked-read sequences')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def qc(inputs, output_dir, min_length, max_length, trim_adapters, deduplicate, deconvolve, extra_params, ignore_bx, threads, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Adapter removal and other FASTQ preprocessing

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    **Standard trimming**
    - a sliding window from front to tail
    - poly-G tail removal

    **Optional quality checks**
    - `-a` remove adapters
      - accepts `auto` to automatically detect adapters or a FASTA file of adapters to remove
    - `-d` finds and removes optical PCR duplicates
      - recommended to skip at this step in favor of barcode-assisted deduplication after alignment
    - `-c` resolves barcodes shared between unrelated sequences
      - off by default, activated with [4 integers](https://github.com/RolandFaure/QuickDeconvolution?tab=readme-ov-file#usage), separated by commas. `21,40,3,0` would be the QuickDeconvolution defaults
      - use `harpy deconvolve` to perform this task separately
    """
    workflow = "qc"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    fqlist, sample_count = parse_fastq_inputs(inputs, "INPUTS")
    if trim_adapters:
        if trim_adapters != "auto":
            if not os.path.exists(trim_adapters):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not exist. Please check the spelling or verify the location of the file.")
            if not os.access(trim_adapters, os.R_OK):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not have read permissions. Please modify the persmissions of the file to grant read access.")
            check_fasta(trim_adapters)
            trim_adapters = filepath(trim_adapters)
    else:
        trim_adapters = False

    ## setup workflow ##
    command, command_rel = setup_snakemake(
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "qc.smk")
    fetch_report(workflowdir, "bx_count.qmd")

    conda_envs = ["qc", "r"]
    configs = {
        "workflow" : workflow,
        "ignore_bx" : ignore_bx,
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
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "reports" : {"skip": skip_reports},
        "inputs" : fqlist
    }
    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:", sample_count),
        ("Trim Adapters:", "yes" if trim_adapters else "no"),
        ("Deduplicate:", "yes" if deduplicate else "no"),
        ("Deconvolve:", "yes" if deconvolve else "no"),
        ("Output Folder:", f"{output_dir}/"),
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/qc.summary")
