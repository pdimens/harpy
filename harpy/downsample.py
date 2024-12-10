"""Downsample a fastq/bam file by barcodes"""

import os
import re
import sys
import yaml
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
            "options": sorted(["--downsample", "--invalid", "--random-seed", "--prefix"]),
        },
        {
            "name": "Workflow Controls",
            "options": ["--quiet", "--snakemake", "--threads", "--help"],
        },
    ]
}
@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/downsample")
@click.option('-d', '--downsample', type = click.IntRange(min = 1), help = 'Downsampling amount')
@click.option('-i', '--invalid', default = 1, show_default = True, type=click.FloatRange(min=0,max=1), help = "Proportion of invalid barcodes to sample")
@click.option('-p', '--prefix', type = click.Path(exists = False), default = "downsampled", show_default = True, help = 'Prefix for output file(s)')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for sampling")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False), nargs=-1)
def downsample(input, prefix, downsample, invalid, random_seed, setup_only, snakemake, threads, quiet):
    """
    Downsample data by barcode

    Downsamples FASTQ or BAM file(s) by barcode in the `BX` tag to keep all reads
    from `-d` barcodes. Input can be:
    - one BAM file
    - two FASTQ files (R1 and R2 of a paired-end read set)
    
    Use `--invalid` to specify the proportion of invalid barcodes to consider for sampling, where `0` will remove all
    invalid barcodes from the barcode pool, `1` will add all invalid barcodes to the barcode pool, and a number in between
    (e.g. `0.25`) will add approximately that proprotion of invalid barcodes into the barcode pool that gets subsampled down
    to `-d` barcodes.
    """
    # validate input files as either 1 bam or 2 fastq
    if len(input) > 2:
        raise click.BadParameter('inputs must be 1 BAM file or 2 FASTQ files.')
    if len(input) == 1:
        if not input[0].lower().endswith(".bam"):
            raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.')
    if len(input) == 2:
        if input[0] == input[1]:
            raise click.BadParameter('the two input files cannot be identical')
        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        for i in input:
            if not re.search(re_ext, i):
                raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.')            
    workflow = "downsample"
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/{workflow}.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if snakemake:
        command += snakemake

    os.makedirs(workflowdir, exist_ok=True)
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    fetch_rule(workflowdir, f"{workflow}.smk")
    sm_log = snakemake_log(output_dir, workflow)
    configs = {
        "workflow": workflow,
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "prefix" :  prefix,
        "downsample" :  downsample,
        "invalid_proportion" : invalid,       
        **({"random_seed" : random_seed} if random_seed else {}),
        "workflow_call" : command.rstrip(),
        "inputs": [i.as_posix() for i in input]
        }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column("value", justify="left")
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Downsample to:", f"{downsample} barcodes")
    start_text.add_row("Invalid Proportion:", f"{invalid}")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, workflow, start_text, output_dir, sm_log, quiet, f"workflow/{workflow}.summary")
