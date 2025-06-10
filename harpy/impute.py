"""Harpy imputation workflow"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._cli_types_generic import convert_to_int, HPCProfile, InputFile, SnakemakeParams
from ._cli_types_params import StitchParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, fetch_report, instantiate_dir, setup_snakemake, write_workflow_config
from ._parsers import parse_alignment_inputs, biallelic_contigs, parse_impute_regions, contigs_from_vcf
from ._printing import workflow_info, print_error, print_solution
from ._validations import vcf_sample_match, check_impute_params, validate_bam_RG

docstring = {
        "harpy impute": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--grid-size", "--region", "--vcf-samples"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/impute/")
@click.option('-x', '--extra-params', type = StitchParams(), help = 'Additional STITCH parameters, in quotes')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Impute", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(2,999, clamp = True), help = 'Number of threads to use')
@click.option('-r', '--region', type = str, help = 'Specific region to impute')
@click.option('-g', '--grid-size', show_default = True, default = 1, type = click.IntRange(min = 1), help = 'Perform imputation in windows of a specific size, instead of per-SNP (default)')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for imputation rather than those found the inputs')
@click.argument('parameters', required = True, type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), nargs=1)
@click.argument('vcf', required = True, type = click.Path(exists=True, readable=True, dir_okay = False, resolve_path=True), nargs=1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def impute(parameters, vcf, inputs, output_dir, region, grid_size, threads, vcf_samples, extra_params, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Impute genotypes using variants and alignments
    
    Provide the parameter file followed by the input VCF and the input alignment files/directories (`.bam`) at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.
    
    Use `harpy template` to generate one and adjust it for your study. Set a `--grid-size` (in base pairs)
    to significantly reduce computation time and memory usage at the cost of minor accuracy loss.
    The `--vcf-samples` option considers only the samples present in your input `VCF` file rather than all
    the samples identified in `INPUTS`. Use `--region` to only impute a specific genomic region, given as
    `contig:start-end-buffer`, otherwise all contigs will be imputed. If providing additional STITCH arguments, they
    must be in quotes and in the `--option=value` format, without spaces (e.g. `"--switchModelIteration=39"`).
    """
    workflow = "impute"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    params = check_impute_params(parameters)
    bamlist, n = parse_alignment_inputs(inputs, "INPUTS")
    validate_bam_RG(bamlist, threads, quiet)
    samplenames = vcf_sample_match(vcf, bamlist, vcf_samples)
    biallelic_file, biallelic_names, n_biallelic = biallelic_contigs(vcf, workflowdir)
    if region:
        contig, start,end, buffer = parse_impute_regions(region, vcf)
        if contig not in biallelic_names:
            print_error(
                "missing contig",
                f"The [bold yellow]{contig}[/] contig given in [blue]{region}[/] is not in the list of contigs identified to have at least 2 biallelic SNPs, therefore it cannot be processed."
            )
            print_solution(f"Restrict the contigs provided to [bold green]--regions[/] to those with at least 2 biallelic SNPs. The contigs Harpy found with at least 2 biallelic can be reviewed in [blue]{biallelic_file}[/].")
            sys.exit(1)
    
    ## setup workflow ##
    command,command_rel = setup_snakemake(
        workflow,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "impute.smk")
    fetch_report(workflowdir, "impute.qmd")
    fetch_report(workflowdir, "stitch_collate.qmd")

    conda_envs = ["r", "stitch"]
    configs = {
        "workflow" : workflow,
        **({'stitch_extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "samples_from_vcf" : vcf_samples,
        **({'region': region} if region else {}),
        "reports" : {"skip": skip_reports},
        "grid_size": grid_size,
        "stitch_parameters" : params,
        "inputs" : {
            "paramfile" : parameters,
            "variantfile" : vcf,
            **({"biallelic_contigs" : biallelic_file} if not region else {}), 
            "alignments" : bamlist
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Input VCF:", os.path.basename(vcf)),
        ("Samples in VCF:", len(samplenames)),
        ("Alignment Files:", n),
        ("Parameter File:", os.path.basename(parameters)),
        ("Contigs:", f"{n_biallelic} [dim](with at least 2 biallelic SNPs)") if not region else ("Target Region:", region),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/impute.summary")
