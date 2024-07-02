"""Harpy sequence adapter trimming and quality control"""

import os
import sys
import subprocess
import rich_click as click

from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_report, fetch_rule
from .fileparsers import parse_fastq_inputs
from .printfunctions import print_onstart

docstring = {
    "harpy qc": [
        {
            "name": "Parameters",
            "options": ["--min-length", "--max-length", "--deduplicate", "--deconvolve", "--deconvolve-params", "--ignore-adapters", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/qc")
@click.option('-n', '--min-length', default = 30, show_default = True, type=int, help = 'Discard reads shorter than this length')
@click.option('-m', '--max-length', default = 150, show_default = True, type=int, help = 'Maximum length to trim sequences down to')
@click.option('-a', '--ignore-adapters', is_flag = True, show_default = False, default = False, help = 'Skip adapter trimming')
@click.option('-d', '--deduplicate', is_flag = True, show_default = True, default = False, help = 'Remove PCR duplicate sequences')
@click.option('-c', '--deconvolve', is_flag = True, show_default = True, default = False, help = 'Resolve barcode clashes between reads from different molecules.')
@click.option('-p', '--deconvolve-params', type = (int,int,int,int), show_default = True, default = (21,40,3,0), help = ' Accepts the QuickDeconvolution parameters for k,w,d,a, in that order')
@click.option('-x', '--extra-params', type = str, help = 'Additional Fastp parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "QC", show_default=True, help = 'Output directory name')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Directory with HPC submission config.yaml file')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--config-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Create the config.yaml file and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def qc(inputs, output_dir, min_length, max_length, ignore_adapters, deduplicate, deconvolve, deconvolve_params, extra_params, threads, snakemake, skipreports, quiet, hpc, conda, config_only):
    """
    Remove adapters and quality-control sequences

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    By default, adapters will be automatically detected and removed (can be disabled with `-a`).
    Use `--deconvolve` to also resolve barcode clashing that may occur by unrelated sequences having the same barcode.
    The parameters are described [here](https://github.com/RolandFaure/QuickDeconvolution?tab=readme-ov-file#usage). You
    can also use the `harpy deconvolve` worklfow to perform this task separately.
    The input reads will be quality trimmed using:
    - a sliding window from front to tail
    - poly-G tail removal
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock  --software-deployment-method {sdm} --conda-prefix .snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/qc.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    os.makedirs(workflowdir, exist_ok=True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    fetch_rule(workflowdir, "qc.smk")
    fetch_report(workflowdir, "BxCount.Rmd")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: qc\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"skip_adapter_trim: {ignore_adapters}\n")
        config.write(f"deduplicate: {deduplicate}\n")
        if deconvolve:
            config.write("deconvolve:\n")
            k,w,d,a = deconvolve_params
            config.write(f"  kmer_length: {k}\n")
            config.write(f"  window_size: {w}\n")
            config.write(f"  density: {d}\n")
            config.write(f"  dropout: {a}\n")
        config.write(f"min_len: {min_length}\n")
        config.write(f"max_len: {max_length}\n")
        config.write(f"extra: {extra_params}\n") if extra_params else None
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        for i in fqlist:
            config.write(f"  - {i}\n")
    if config_only:
        sys.exit(0)

    print_onstart(
        f"Samples: {sample_count}\nOutput Directory: {output_dir}/",
        "qc"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)
