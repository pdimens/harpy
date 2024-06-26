"""Harpy sequence adapter trimming and quality control"""

import os
import sys
import subprocess
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_report, fetch_rule, fetch_script
from .fileparsers import parse_fastq_inputs
from .printfunctions import print_onstart

docstring = {
    "harpy qc": [
        {
            "name": "Parameters",
            "options": ["--max-length", "--ignore-adapters", "--extra-params"],
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
@click.option('-x', '--extra-params', type = str, help = 'Additional Fastp parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "QC", show_default=True, help = 'Name of output directory')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def qc(inputs, output_dir, min_length, max_length, ignore_adapters, extra_params, threads, snakemake, skipreports, quiet, hpc, conda, print_only):
    """
    Remove adapters and quality trim sequences

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    By default, adapters will be automatically detected and removed (can be disabled with `-a`).
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

    if print_only:
        click.echo(command)
        sys.exit(0)

    os.makedirs(workflowdir, exist_ok=True)
    fqlist, sample_count = parse_fastq_inputs(inputs)
    fetch_rule(workflowdir, "qc.smk")
    fetch_report(workflowdir, "BxCount.Rmd")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: qc\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"skip_adapter_trim: {ignore_adapters}\n")
        config.write(f"min_len: {min_length}\n")        
        config.write(f"max_len: {max_length}\n")
        config.write(f"extra: {extra_params}\n") if extra_params else None
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        for i in fqlist:
            config.write(f"  - {i}\n")

    print_onstart(
        f"Samples: {sample_count}\nOutput Directory: {output_dir}/",
        "qc"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)
