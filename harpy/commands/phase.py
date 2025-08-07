"""Harpy haplotype phasing workflow"""

import os
import rich_click as click
from harpy.common.cli_types_generic import ContigList, HPCProfile, InputFile, SnakemakeParams
from harpy.common.cli_types_params import HapCutParams
from harpy.common.misc import container_ok
from harpy.common.parsers import parse_alignment_inputs
from harpy.common.printing import workflow_info
from harpy.common.validations import check_fasta, vcf_sample_match, validate_bam_RG, vcf_contig_match
from harpy.common.workflow import Workflow

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
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = click.IntRange(min = 100), help = 'Distance cutoff to split molecules (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Phase", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prune-threshold', default = 30, show_default = True, type = click.IntRange(0,100, clamp = True), help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(2, 999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
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
    workflow = Workflow("phase", "phase.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["hapcut.qmd"]
    workflow.conda = ["phase", "r"]

    ## checks and validations ##
    bamlist, n = parse_alignment_inputs(inputs, "INPUTS")
    samplenames = vcf_sample_match(vcf, bamlist, vcf_samples)
    validate_bam_RG(bamlist, threads, quiet)
    if reference:
        check_fasta(reference)
    if contigs:
        vcf_contig_match(contigs, vcf)

    workflow.config = {
        "workflow" : workflow.name,
        "prune" : prune_threshold,
        "samples_from_vcf" : vcf_samples,
        "barcodes": {
            "ignore" : ignore_bx,
            "distance_threshold" : molecule_distance,
        },
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
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

    workflow.start_text = workflow_info(
        ("Input VCF:", os.path.basename(vcf)),
        ("Samples in VCF:", len(samplenames)),
        ("Alignment Files:", n),
        ("Phase Indels:", "yes" if reference else "no"),
        ("Reference:", os.path.basename(reference)) if reference else None,
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)
