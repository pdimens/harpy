"""Downsample a fastq/bam file by barcodes"""

import os
import re
import sys
import yaml
import shutil
from pathlib import Path
import rich_click as click
from ._cli_types_generic import convert_to_int, SnakemakeParams
from ._launch import launch_snakemake, SNAKEMAKE_CMD
from ._misc import fetch_rule, snakemake_log
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
@click.option('-d', '--downsample', type = click.IntRange(min = 1), help = 'Downsampling amount')
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
    Downsample reads by barcode from FASTQ or BAM files.
    
    This function down-samples input read files based on the terminal "BX:Z" barcode tag
    to retain all reads associated with a specified number of barcodes. Input must be a single
    BAM file or a pair of FASTQ files (for paired-end reads). If the "BX:Z" tag is not terminal,
    use the bx_to_end.py utility to reposition it.
    
    The invalid barcode proportion is configurable: a value of 0 excludes all invalid barcodes,
    1 includes all invalid barcodes, and a fractional value between 0 and 1 retains that fraction.
    An optional HPC submission YAML configuration file may be provided via the hpc parameter.
    When specified, it is copied into the workflow directory and applied as a workflow profile.
    If setup_only is enabled, the workflow environment is configured and the function exits
    without executing the workflow.
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
    command = f'{SNAKEMAKE_CMD} --cores {threads}'
    command += f" --snakefile {workflowdir}/{workflow}.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

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
        "inputs": [Path(i).resolve().as_posix() for i in input]
        }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Output Folder:", output_dir + "/"),
        ("Downsample to:", f"{downsample} barcodes"),
        ("Invalid Proportion:", invalid),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, workflow, start_text, output_dir, sm_log, quiet, f"workflow/{workflow}.summary")
