"""Separate barcodes by unique molecule"""

import os
import sys
import yaml
import rich_click as click
from ._cli_types_generic import convert_to_int, HPCProfile, SnakemakeParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake, SNAKEMAKE_CMD
from ._misc import fetch_rule, snakemake_log
from ._parsers import parse_fastq_inputs
from ._printing import workflow_info

docstring = {
    "harpy deconvolve": [
        {
            "name": "Parameters",
            "options": ["--density", "--dropout", "--kmer-length", "--window-size"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/qc")
@click.option('-k', '--kmer-length', default = 21, show_default = True, type=int, help = 'Size of kmers')
@click.option('-w', '--window-size', default = 40, show_default = True, type=int, help = 'Size of window guaranteed to contain at least one kmer')
@click.option('-d', '--density', default = 3, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'On average, 1/2^d kmers are indexed')
@click.option('-a', '--dropout', default = 0, show_default = True, type = click.IntRange(min = 0, max_open = True), help = 'Minimum cloud size to deconvolve')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max = 999), help = 'Number of threads to use')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Deconvolve", show_default=True,  help = 'Output directory name')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def deconvolve(inputs, output_dir, kmer_length, window_size, density, dropout, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Resolve shared barcodes between unrelated molecules

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    The term "cloud" refers to the collection of all sequences that feature the same barcode. By default,
    `dropout` is set to `0`, meaning it will consider all barcodes, even clouds with singleton.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/deconvolve.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        command += f" --workflow-profile {hpc}"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(workflowdir, exist_ok=True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    fetch_rule(workflowdir, "deconvolve.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "deconvolve")
    conda_envs = ["qc"]
    configs = {
        "workflow": "deconvolve",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "kmer_length" : kmer_length,       
        "window_size" : window_size,
        "density" :  density,
        "dropout" :  dropout,
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "inputs": [i.as_posix() for i in fqlist]
        }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:", sample_count),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "deconvolve", start_text, output_dir, sm_log, quiet, "workflow/deconvolve.summary")
