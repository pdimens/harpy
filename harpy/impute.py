"""Harpy imputation workflow"""

import os
import sys
import yaml
from pathlib import Path
import rich_click as click
from ._cli_types_generic import convert_to_int, HPCProfile, InputFile, SnakemakeParams
from ._cli_types_params import StitchParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake, SNAKEMAKE_CMD
from ._misc import fetch_rule, fetch_report, fetch_script, snakemake_log
from ._parsers import parse_alignment_inputs, biallelic_contigs
from ._printing import workflow_info
from ._validations import vcf_sample_match, check_impute_params, validate_bam_RG

docstring = {
        "harpy impute": [
        {
            "name": "Parameters",
            "options": ["--extra-params",  "--parameters", "--vcf", "--vcf-samples"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/impute/")
@click.option('-x', '--extra-params', type = StitchParams(), help = 'Additional STITCH parameters, in quotes')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Impute", show_default=True,  help = 'Output directory name')
@click.option('-p', '--parameters', required = True, type=click.Path(exists=True, dir_okay=False, readable=True), help = 'STITCH parameter file (tab-delimited)')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-v', '--vcf', required = True, type = InputFile("vcf", gzip_ok = False), help = 'Path to BCF/VCF file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = 'Verbosity of output. `0` shows all output, `1` shows single progress bar, `2` suppressess all output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for imputation rather than those found the inputs')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def impute(inputs, output_dir, parameters, threads, vcf, vcf_samples, extra_params, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Impute genotypes using variants and alignments
    
    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.
    
    Requires a parameter file, use **harpy stitchparams** to generate one and adjust it for your study.
    The `--vcf-samples` option considers only the samples present in your input `--vcf` file rather than all
    the samples identified in `INPUT`. If providing additional STITCH arguments, they must be in quotes and 
    in the `--option=value` format, without spaces between `--option` and `value` (e.g. `"--switchModelIteration=39"`).
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/impute.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        command += f" --workflow-profile {hpc}"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    params = check_impute_params(parameters)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_bam_RG(bamlist, threads, quiet)
    samplenames = vcf_sample_match(vcf, bamlist, vcf_samples)
    biallelic, n_biallelic = biallelic_contigs(vcf, f"{workflowdir}")

    fetch_rule(workflowdir, "impute.smk")
    fetch_report(workflowdir, "impute.qmd")
    fetch_report(workflowdir, "stitch_collate.qmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "impute")
    conda_envs = ["r", "stitch"]
    configs = {
        "workflow" : "impute",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "samples_from_vcf" : vcf_samples,
        **({'stitch_extra': extra_params} if extra_params else {}),
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {"skip": skip_reports},
        "stitch_parameters" : params,
        "inputs" : {
            "paramfile" : Path(parameters).resolve().as_posix(),
            "variantfile" : Path(vcf).resolve().as_posix(),
            "biallelic_contigs" : Path(biallelic).resolve().as_posix(),
            "alignments" : [i.as_posix() for i in bamlist]
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))
    
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Input VCF:", vcf),
        ("VCF Samples:", len(samplenames)),
        ("Alignment Files:", n),
        ("Parameter File:", parameters),
        ("Contigs:", f"{n_biallelic} [dim](with at least 2 biallelic SNPs)"),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "impute", start_text, output_dir, sm_log, quiet, "workflow/impute.summary")
