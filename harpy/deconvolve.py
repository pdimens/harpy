"""Separate barcodes by unique molecule"""

import os
import sys
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, snakemake_log, launch_snakemake
from .fileparsers import parse_fastq_inputs

docstring = {
    "harpy deconvolve": [
        {
            "name": "Parameters",
            "options": ["--density", "--dropout", "--kmer-length", "--window-size"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/qc")
@click.option('-k', '--kmer-length', default = 21, show_default = True, type=int, help = 'Size of kmers')
@click.option('-w', '--window-size', default = 40, show_default = True, type=int, help = 'Size of window guaranteed to contain at least one kmer')
@click.option('-d', '--density', default = 3, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'On average, 1/2^d kmers are indexed')
@click.option('-a', '--dropout', default = 0, show_default = True, type = click.IntRange(min = 0, max_open = True), help = 'Minimum cloud size to deconvolve')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Deconvolve", show_default=True,  help = 'Output directory name')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--config-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Create the config.yaml file and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def deconvolve(inputs, output_dir, kmer_length, window_size, density, dropout, threads, snakemake, quiet, hpc, conda, config_only):
    """
    Resolve clashing barcodes from different molecules

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    The term "cloud" refers to the collection of all sequences that feature the same barcode. By default,
    `dropout` is set to `0`, meaning it will consider all barcodes, even clouds with singleton.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock  --software-deployment-method {sdm} --conda-prefix .snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/deconvolve.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake is not None:
        command += snakemake

    os.makedirs(workflowdir, exist_ok=True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    fetch_rule(workflowdir, "deconvolve.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "deconvolve")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: deconvolve\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"kmer_length: {kmer_length}\n")       
        config.write(f"window_size: {window_size}\n")
        config.write(f"density: {density}\n")
        config.write(f"dropout: {dropout}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        for i in fqlist:
            config.write(f"  - {i}\n")
    if config_only:
        sys.exit(0)

    generate_conda_deps()
    start_text =  f"Samples: {sample_count}\nOutput Directory: {output_dir}/\nLog: {sm_log}"
    launch_snakemake(command, "deconvolve", start_text, output_dir, sm_log, quiet)