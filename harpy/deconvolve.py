"""Separate barcodes by unique molecule"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._cli_types_generic import HPCProfile, SnakemakeParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, instantiate_dir, setup_snakemake, write_workflow_config
from ._parsers import parse_fastq_inputs
from ._printing import workflow_info

docstring = {
    "harpy deconvolve": [
        {
            "name": "Parameters",
            "options": ["--density", "--dropout", "--kmer-length", "--window-size"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/qc")
@click.option('-k', '--kmer-length', default = 21, show_default = True, type=int, help = 'Size of kmers')
@click.option('-w', '--window-size', default = 40, show_default = True, type=int, help = 'Size of window guaranteed to contain at least one kmer')
@click.option('-d', '--density', default = 3, show_default = True, type = click.IntRange(min = 1), help = 'On average, 1/2^d kmers are indexed')
@click.option('-a', '--dropout', default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Minimum cloud size to deconvolve')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1,999, clamp = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Deconvolve", show_default=True,  help = 'Output directory name')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def deconvolve(inputs, output_dir, kmer_length, window_size, density, dropout, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Resolve barcode sharing in unrelated molecules

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    The term "cloud" refers to the collection of all sequences that feature the same barcode. By default,
    `dropout` is set to `0`, meaning it will consider all barcodes, even clouds with singleton.
    """
    workflow = "deconvolve"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    fqlist, sample_count = parse_fastq_inputs(inputs, "INPUTS")
    
    ## setup workflow ##
    command,command_rel = setup_snakemake(
        workflow,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "deconvolve.smk")

    conda_envs = ["qc"]
    configs = {
        "workflow": workflow,
        "kmer_length" : kmer_length,       
        "window_size" : window_size,
        "density" :  density,
        "dropout" :  dropout,
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "inputs": fqlist
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:", sample_count),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/deconvolve.summary")
