"""Harpy imputation workflow"""

import os
import sys
import yaml
import shutil
from pathlib import Path
import rich_click as click
from ._cli_types_generic import convert_to_int, HPCProfile, InputFile, SnakemakeParams
from ._cli_types_params import StitchParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, fetch_report, fetch_script, snakemake_log, write_snakemake_config, write_workflow_config
from ._parsers import parse_alignment_inputs, biallelic_contigs, parse_impute_regions, contigs_from_vcf
from ._printing import workflow_info, print_error, print_solution
from ._validations import vcf_sample_match, check_impute_params, validate_bam_RG

docstring = {
        "harpy impute": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--regions", "--vcf-samples"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/impute/")
@click.option('-x', '--extra-params', type = StitchParams(), help = 'Additional STITCH parameters, in quotes')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Impute", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-r', '--regions', type = click.Path(exists=True, dir_okay=False, readable=True), help = 'BED file of regions to impute')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for imputation rather than those found the inputs')
@click.argument('parameters', required = True, type=click.Path(exists=True, dir_okay=False, readable=True), nargs=1)
@click.argument('vcf', required = True, type = click.Path(exists=True, readable=True, dir_okay = False), nargs=1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def impute(parameters, vcf, inputs, output_dir, regions, threads, vcf_samples, extra_params, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Impute genotypes using variants and alignments
    
    Provide the input VCF followed by the input alignment files (`.bam`) and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.
    
    Requires a parameter file, use **harpy template** to generate one and adjust it for your study.
    The `--vcf-samples` option considers only the samples present in your input `VCF` file rather than all
    the samples identified in `INPUTS`. Use `--regions` to only impute specific genomic regions, given as a BED file
    of whitespace-delimited `chromosome start end`. If providing additional STITCH arguments, they must be in quotes and 
    in the `--option=value` format, without spaces between `--option` and `value` (e.g. `"--switchModelIteration=39"`).
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    os.makedirs(workflowdir, exist_ok= True)

    ## checks and validations ##
    params = check_impute_params(parameters)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_bam_RG(bamlist, threads, quiet)
    samplenames = vcf_sample_match(vcf, bamlist, vcf_samples)
    biallelic_file, biallelic_names, n_biallelic = biallelic_contigs(vcf, workflowdir)
    if regions:
        target_regions = parse_impute_regions(regions, vcf)
        for i in target_regions:
            # why did I write this? What is it protecting us from?
            #chrm_name = i.split("_")[0]
            if i not in biallelic_names:
                print_error(
                    "missing contig",
                    f"The [bold yellow]{i}[/bold yellow] contig given in [blue]{regions}[/blue] is not in the list of contigs identified to have at least 2 biallelic SNPs, therefore it cannot be processed."
                )
                print_solution(f"Restrict the contigs provided to [bold green]--regions[/bold green] to those with at least 2 biallelic SNPs. The contigs Harpy found with at least 2 biallelic can be reviewed in [blue]{biallelic_file}[/blue].")
                sys.exit(1)
    else:
        # get the contigs and their lengths from the VCF file
        target_regions = contigs_from_vcf(vcf)
        # filter for only biallelic contigs and prepend with 1 as start position
        target_regions = {_contig: [1,target_regions[_contig]] for _contig in biallelic_contigs}
    with open(f"{workflowdir}/regions.bed", "w") as workflow_bed:
        for k,v in target_regions:
            workflow_bed.write(f"{k}\t{v[0]}\t{v[1]}\n")

    ## setup workflow ##
    write_snakemake_config("conda" if not container else "conda apptainer", output_dir)
    command = f"snakemake --cores {threads} --snakefile {workflowdir}/impute.smk"
    command += f" --configfile {workflowdir}/config.harpy.yaml --profile {workflowdir}"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    fetch_rule(workflowdir, "impute.smk")
    fetch_report(workflowdir, "impute.qmd")
    fetch_report(workflowdir, "stitch_collate.qmd")

    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "impute")
    conda_envs = ["r", "stitch"]
    configs = {
        "workflow" : "impute",
        "snakemake_log" : sm_log,
        "samples_from_vcf" : vcf_samples,
        **({'stitch_extra': extra_params} if extra_params else {}),
        "snakemake_command" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {"skip": skip_reports},
        "stitch_parameters" : params,
        "regions" : target_regions,
        "inputs" : {
            "paramfile" : Path(parameters).resolve().as_posix(),
            "variantfile" : Path(vcf).resolve().as_posix(),
            "biallelic_contigs" : Path(biallelic_file).resolve().as_posix(),
            "alignments" : [i.as_posix() for i in bamlist]
        }
    }

    write_workflow_config(config, workflowdir)
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
