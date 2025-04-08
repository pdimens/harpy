"""Harpy sequence adapter trimming and quality control"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_report, fetch_rule, snakemake_log, write_snakemake_config, write_workflow_config
from ._cli_types_generic import convert_to_int, HPCProfile, SnakemakeParams
from ._cli_types_params import FastpParams
from ._misc import filepath
from ._parsers import parse_fastq_inputs
from ._printing import workflow_info
from ._validations import check_fasta

docstring = {
    "harpy qc": [
        {
            "name": "Parameters",
            "options": ["--deconvolve", "--deconvolve-params", "--deduplicate", "--extra-params", "--min-length", "--max-length", "--trim-adapters"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--ignore-bx", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/qc")
@click.option('-c', '--deconvolve', type = click.IntRange(0,100, clamp = True), nargs=4, help = "`k` `w` `d` `a` QuickDeconvolution parameters")
@click.option('-d', '--deduplicate', is_flag = True, default = False, help = 'Identify and remove PCR duplicates')
@click.option('-x', '--extra-params', type = FastpParams(), help = 'Additional Fastp parameters, in quotes')
@click.option('-m', '--max-length', default = 150, show_default = True, type=int, help = 'Maximum length to trim sequences down to')
@click.option('-n', '--min-length', default = 30, show_default = True, type=int, help = 'Discard reads shorter than this length')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "QC", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-a', '--trim-adapters', type = str, help = 'Detect and trim adapters')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--ignore-bx',  is_flag = True, default = False, help = 'Ignore parts of the workflow specific to linked-read sequences')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def qc(inputs, output_dir, min_length, max_length, trim_adapters, deduplicate, deconvolve, extra_params, ignore_bx, threads, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Remove adapters and quality-control sequences

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    **Standard trimming**
    - a sliding window from front to tail
    - poly-G tail removal

    **Optional quality checks**
    - `-a` remove adapters
      - accepts `auto` to automatically detect adapters or a FASTA file of adapters to remove
    - `-d` finds and removes optical PCR duplicates
      - recommended to skip at this step in favor of barcode-assisted deduplication after alignment
    - `-c` resolves barcodes shared between unrelated sequences
      - off by default, activated with [4 integers](https://github.com/RolandFaure/QuickDeconvolution?tab=readme-ov-file#usage), separated by spaces. `21 40 3 0` would be the QuickDeconvolution defaults
      - use `harpy deconvolve` to perform this task separately
    """
    ## checks and validations ##
    fqlist, sample_count = parse_fastq_inputs(inputs)
    if trim_adapters:
        if trim_adapters != "auto":
            if not os.path.exists(trim_adapters):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not exist. Please check the spelling or verify the location of the file.")
            if not os.access(trim_adapters, os.R_OK):
                raise click.BadParameter(f"--trim-adapters was given {trim_adapters}, but that file does not have read permissions. Please modify the persmissions of the file to grant read access.")
            check_fasta(trim_adapters)
            trim_adapters = filepath(trim_adapters)
    else:
        trim_adapters = False

    ## setup workflow ##
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    write_snakemake_config("conda" if not container else "conda apptainer", output_dir)
    command = f"snakemake --cores {threads} --snakefile {workflowdir}/qc.smk"
    command += f" --configfile {workflowdir}/config.harpy.yaml --profile {workflowdir}"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    fetch_rule(workflowdir, "qc.smk")
    fetch_report(workflowdir, "bx_count.qmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "qc")
    conda_envs = ["qc", "r"]
    configs = {
        "workflow" : "qc",
        "snakemake_log" : sm_log,
        "ignore_bx" : ignore_bx,
        "trim_adapters" : trim_adapters,
        "deduplicate" : deduplicate,
        "min_len" : min_length,
        "max_len" : max_length,
        **({'extra': extra_params} if extra_params else {}),
        **({'deconvolve': {
            "kmer_length" : deconvolve[0],
            "window_size" : deconvolve[1],
            "density" : deconvolve[2],
            "dropout" : deconvolve[3]
        }} if deconvolve else {}),
        "snakemake_command" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {"skip": skip_reports},
        "inputs" : fqlist
    }
    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:", sample_count),
        ("Trim Adapters:", "yes" if trim_adapters else "no"),
        ("Deduplicate:", "yes" if deduplicate else "no"),
        ("Deconvolve:", "yes" if deconvolve else "no"),
        ("Output Folder:", f"{output_dir}/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "qc", start_text, output_dir, sm_log, quiet, "workflow/qc.summary")
