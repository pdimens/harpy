"""Harpy haplotype phasing workflow"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, fetch_report, instantiate_dir, setup_snakemake, write_workflow_config
from ._cli_types_generic import convert_to_int, ContigList, HPCProfile, InputFile, SnakemakeParams
from ._cli_types_params import HapCutParams
from ._parsers import parse_alignment_inputs
from ._printing import workflow_info
from ._validations import check_fasta, vcf_sample_match, validate_bam_RG, vcf_contig_match

docstring = {
        "harpy phase": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--reference", "--ignore-bx", "--molecule-distance", "--prune-threshold", "--vcf-samples"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        }
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/phase")
@click.option('-x', '--extra-params', type = HapCutParams(), help = 'Additional HapCut2 parameters, in quotes')
@click.option('-r', '--reference', type=InputFile("fasta", gzip_ok = True), help = 'Path to reference genome if wanting to also extract reads spanning indels')
@click.option('-b', '--ignore-bx',  is_flag = True, show_default = True, default = False, help = 'Ignore barcodes when phasing')
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Distance cutoff to split molecules (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Phase", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prune-threshold', default = 7, show_default = True, type = click.IntRange(0,100, clamp = True), help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(2, 999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for phasing rather than those found the inputs')
@click.argument('vcf', required = True, type = InputFile("vcf", gzip_ok = False), nargs = 1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def phase(vcf, inputs, output_dir, threads, molecule_distance, prune_threshold, vcf_samples, reference, snakemake, extra_params, ignore_bx, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Phase SNPs into haplotypes

    Provide the vcf file followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/myotis*.bam`), or both.
    
    You may choose to omit barcode information with `--ignore-bx`, although it's usually
    better to include that information. Use `--vcf-samples` to phase only
    the samples present in your input `VCF` file rather than all the samples present in
    the `INPUT` alignments.
    """
    workflow = "phase"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    bamlist, n = parse_alignment_inputs(inputs)
    samplenames = vcf_sample_match(vcf, bamlist, vcf_samples)
    validate_bam_RG(bamlist, threads, quiet)
    if reference:
        check_fasta(reference)
    if contigs:
        vcf_contig_match(contigs, vcf)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        workflow,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "phase.smk")
    fetch_report(workflowdir, "hapcut.qmd")

    conda_envs = ["phase", "r"]
    configs = {
        "workflow" : workflow,
        "prune" : prune_threshold/100,
        "samples_from_vcf" : vcf_samples,
        "barcodes": {
            "ignore" : ignore_bx,
            "distance_threshold" : molecule_distance,
        },
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "variantfile" : vcf,
            **({'reference': reference} if reference else {}),
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
        ("Phase Indels:", "yes" if reference else "no"),
        ("Reference:", os.path.basename(reference)) if reference else None,
        ("Output Folder:", os.path.basename(output_dir) + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/phase.summary")
