"""Downsample a fastq/bam file by barcodes"""

import os
import re
import sys
import yaml
from pathlib import Path
from rich import box
from rich.table import Table
import rich_click as click
from ._cli_types_generic import SnakemakeParams
from ._launch import launch_snakemake
from ._misc import fetch_rule, snakemake_log

docstring = {
    "harpy downsample": [
        {
            "name": "Parameters",
            "options": sorted(["--downsample", "--invalid", "--bx-tag", "--random-seed"]),
        },
        {
            "name": "Workflow Controls",
            "options": ["--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/downsample")
@click.option('-d', '--downsample', type = click.FloatRange(min = 0.000001), help = 'Downsampling amount')
@click.option('-i', '--invalid', default = "keep", show_default = True, type=click.Choice( ["keep","drop", "downsample"]), help = "Strategy to handle invalid/missing barcodes")
@click.option('-b', '--bx-tag', type = str, default = "BX", show_default=True, help = "The header tag with the barcode")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Downsample", show_default=True,  help = 'Output directory name')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for sampling")
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('input', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def downsample(input, output_dir, downsample, invalid, bx_tag, random_seed, threads, snakemake, quiet, setup_only):
    """
    Downsample data by barcode

    Downsamples FASTQ or BAM file(s) by barcode using one of the
    two methods described below. Input files must be either:
    - one BAM file
    - two FASTQ files (R1 and R2 of a paired-end read set)
    
    ### Downsampling Methods
    | -d value | unpaired | function |
    |:---:|:---|:---|
    | < 1 | dropped | keep `-d` proportion of reads for every barcode |
    | > 1 | kept | keep all reads from `-d` barcodes |

    Use `--invalid` to specify how to handle invalid barcodes:
    - `keep`: keep all invalid/missing barcodes
    - `drop`: don't output any invalid/missing barcodes
    - `downsample`: downsample invalid/missing barcodes
      - (`-d < 1` only)
    """
    # validate input files as either 1 bam or 2 fastq
    if len(bx_tag) != 2:
        raise click.BadParameter(f'\'{bx_tag}\' is not a valid SAM tag. Tags for --bx-tag must be exactly 2 characters, e.g. "BX"')
    if len(input) > 2:
        raise click.BadParameter('inputs must be 1 BAM file or 2 FASTQ files.')
    if len(input) == 1:
        if not input[0].lower().endswith(".bam"):
            raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.')
    if len(input) == 2:
        if input[0] == input[1]:
            raise click.BadParameter('the two input files cannot be identical')
        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        badfastq = []
        for i in input:
            if not re.search(re_ext, i):
                badfastq.append(i)
        if badfastq:
            raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.')
    
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/downsample.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if snakemake:
        command += snakemake

    os.makedirs(workflowdir, exist_ok=True)
    fetch_rule(workflowdir, "downsample.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "downsample")
    conda_envs = ["qc"]
    configs = {
        "workflow": "downsample",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "downsample" : downsample if downsample < 1 else int(downsample),       
        "invalid" : invalid,
        "bx_tag" :  bx_tag.upper(),
        **({"random_seed" : random_seed} if random_seed else {}),
        "workflow_call" : command.rstrip(),
        "inputs": [Path(i).as_posix() for i in input]
        }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column("value", justify="left")
    start_text.add_row("Downsample:", f"within {bx_tag} tags" if downsample < 1 else f"across {bx_tag} tags")
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, "deconvolve", start_text, output_dir, sm_log, quiet, "workflow/deconvolve.summary")
