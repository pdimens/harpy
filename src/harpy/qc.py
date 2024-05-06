from .helperfunctions import generate_conda_deps, fetch_report, fetch_rule, fetch_script
from .fileparsers import get_samples_from_fastq, parse_fastq_inputs
from .printfunctions import print_onstart
import rich_click as click
import re
import os
import sys
import glob

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/qc")
@click.option('-n', '--min-length', default = 30, show_default = True, type=int, help = 'Discard reads shorter than this length')
@click.option('-m', '--max-length', default = 150, show_default = True, type=int, help = 'Maximum length to trim sequences down to')
@click.option('-a', '--ignore-adapters', is_flag = True, show_default = False, default = False, help = 'Skip adapter trimming')
@click.option('-x', '--extra-params', type = str, help = 'Additional Fastp parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "QC", show_default=True, help = 'Name of output directory')
@click.option('--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate any HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def qc(input, output_dir, min_length, max_length, ignore_adapters, extra_params, threads, snakemake, skipreports, quiet, print_only):
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
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock  --software-deployment-method conda apptainer apptainer --conda-prefix .snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/qc.smk')
    command.append('--configfile')
    command.append(f'{workflowdir}/config.yml')
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit(0)

    _ = parse_fastq_inputs(input, f"{workflowdir}/input")
    sn = get_samples_from_fastq(f"{workflowdir}/input")

    fetch_script(workflowdir, "countBX.py")    
    fetch_rule(workflowdir, "qc.smk")
    fetch_report(workflowdir, "BxCount.Rmd")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"adapters: {ignore_adapters}\n")
        config.write(f"min_len: {min_length}\n")        
        config.write(f"max_len: {max_length}\n")
        config.write(f"extra: {extra_params}\n") if extra_params else None
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Samples: {len(sn)}\nOutput Directory: {output_dir}/",
        "qc"
    )
    return command