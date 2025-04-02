"""Downsample a fastq/bam file by barcodes"""

import os
import re
import sys
import yaml
import shutil
from pathlib import Path
import rich_click as click
from ._cli_types_generic import convert_to_int, SnakemakeParams, HPCProfile
from ._launch import launch_snakemake
from ._misc import fetch_rule, snakemake_log, write_snakemake_config, write_workflow_config
from ._printing import workflow_info

docstring = {
    "harpy downsample": [
        {
            "name": "Parameters",
            "options": sorted(["--downsample", "--invalid", "--random-seed", "--prefix"]),
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}
@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/downsample")
@click.option('-d', '--downsample', type = click.IntRange(min = 1), help = 'Number of barcodes to retain')
@click.option('-i', '--invalid', default = 1, show_default = True, type=click.FloatRange(min=0,max=1), help = "Proportion of invalid barcodes to sample")
@click.option('-p', '--prefix', type = click.Path(exists = False), default = "downsampled", show_default = True, help = 'Prefix for output file(s)')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for sampling")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Downsample", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1,999, clamp = True), help = 'Number of threads to use')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False), nargs=-1)
def downsample(input, invalid, output_dir, prefix, downsample, random_seed, hpc, setup_only, snakemake, threads, quiet):
    """
    Downsample data by barcode
    
    Downsamples FASTQ or BAM file(s) by barcode in the `BX` tag to keep all reads
    from `-d` barcodes. The `BX:Z` tag must be the **last tag** in the FASTQ file(s).
    If the `BX` tag isn't terminal, use `bx_to_end.py` (provided by Harpy) to move
    the tag to the end. Input can be:
    - one BAM file
    - two FASTQ files (R1 and R2 of a paired-end read set)
    
    Use `--invalid` to specify the proportion of invalid barcodes to consider for sampling:
    - `0` removes all invalid barcodes from the barcode pool
    - `1` adds all invalid barcodes to the barcode pool
    - 0<`i`<1 (e.g. `0.25`) keeps that proprotion of invalid barcodes in the barcode pool
    """
    ## checks and validations ##
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

    ## setup workflow ##
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    write_snakemake_config("conda" if not container else "conda apptainer", output_dir)
    command = f"snakemake --cores {threads} --snakefile {workflowdir}/downsample.smk"
    command += f" --configfile {workflowdir}/config.harpy.yaml --profile {workflowdir}"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)

    fetch_rule(workflowdir, f"{workflow}.smk")
    sm_log = snakemake_log(output_dir, workflow)

    configs = {
        "workflow": "downsample",
        "snakemake_log" : sm_log,
        "prefix" :  prefix,
        "downsample" :  downsample,
        "invalid_proportion" : invalid,       
        **({"random_seed" : random_seed} if random_seed else {}),
        "snakemake_command" : command.rstrip(),
        "inputs": [Path(i).resolve().as_posix() for i in input]
    }

    write_workflow_config(conrigs, workflowdir)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Output Folder:", output_dir + "/"),
        ("Downsample to:", f"{downsample} barcodes"),
        ("Invalid Proportion:", invalid),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "downsample", start_text, output_dir, sm_log, quiet, f"workflow/downsample.summary")
