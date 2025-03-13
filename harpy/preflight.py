"""Harpy preflight-check workflows for FASTQ and BAM files"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._cli_types_generic import convert_to_int, HPCProfile, SnakemakeParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake, SNAKEMAKE_CMD
from ._misc import fetch_rule, fetch_report, snakemake_log
from ._parsers import parse_alignment_inputs, parse_fastq_inputs
from ._printing import workflow_info

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def preflight():
    """
    File format checks for haplotag data

    This is useful to make sure your input files are formatted correctly for the processing pipeline 
    before you are surprised by errors hours into an analysis. Provide an additional command `fastq`
    or `bam` to see more information and options.
    """

docstring = {
    "harpy preflight bam": [
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy preflight fastq": [
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preflight/")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Preflight/bam", show_default=True,  help = 'Output directory name')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def bam(inputs, output_dir, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Run validity checks on haplotagged BAM files

    Provide the input alignment (`.bam`) files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/betula*.bam`), or both.
    
    It will check that alignments have BX:Z: tags, that haplotag
    barcodes are properly formatted (`AxxCxxBxxDxx`) and that the filename matches the `@RG ID` tag.
    This **will not** fix your data, but it will report the number of records that feature errors  to help
    you diagnose if file formatting will cause downstream issues. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/preflight_bam.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    fetch_rule(workflowdir, "preflight_bam.smk")
    fetch_report(workflowdir, "preflight_bam.qmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "preflight_bam")
    conda_envs = ["r"]
    configs = {
        "workflow" : "preflight bam",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "inputs" : [i.as_posix() for i in bamlist]
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Alignment Files:", n),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "preflight_bam", start_text, output_dir, sm_log, quiet, "workflow/preflight.bam.summary")

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preflight/")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Preflight/fastq", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def fastq(inputs, output_dir, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Run validity checks on haplotagged FASTQ files.

    Provide the input fastq files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/wombat*.fastq.gz`), or both.
    
    It will check that fastq reads have `BX:Z:` tags, that haplotag
    barcodes are propery formatted (`AxxCxxBxxDxx`) and that the comments in the
    read headers conform to the SAM specification of `TAG:TYPE:VALUE`. This **will not**
    fix your data, but it will report the number of reads that feature errors to help
    you diagnose if file formatting will cause downstream issues. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/preflight_fastq.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fqlist, n = parse_fastq_inputs(inputs)
    fetch_rule(workflowdir, "preflight_fastq.smk")
    fetch_report(workflowdir, "preflight_fastq.qmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "preflight_fastq")
    conda_envs = ["r"]
    configs = {
        "workflow" : "preflight fastq",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "inputs" : [i.as_posix() for i in fqlist]
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("FASTQ Files:", n),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "preflight_fastq", start_text, output_dir, sm_log, quiet, "workflow/preflight.fastq.summary")

preflight.add_command(bam)
preflight.add_command(fastq)
