"""Harpy haplotype phasing workflow"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, SAMfile, VCFfile, FASTAfile
from harpy.validation.fasta import FASTA
from harpy.validation.xam import XAM
from harpy.validation.vcf import VCF
from harpy.common.cli_types_generic import ContigList, SnakemakeParams
from harpy.common.cli_types_params import HapCutParams
from harpy.common.printing import workflow_info
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow

@click.group(context_settings={"help_option_names" : ['--help']})
def phase():
    """
    Phase SNPs or alignments

    Provide an additional subcommand `snp` or `bam` for more information the phasing workflows.
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/phase")
@click.option('-x', '--extra-params', panel = "Parameters", type = HapCutParams(), help = 'Additional whatshap haplotag parameters, in quotes')
@click.option('-d', '--molecule-distance', panel = "Parameters", default = 100000, show_default = True, type = click.IntRange(min = 100), help = 'Distance cutoff to split molecules (bp)')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Phase/bam", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(3, 999, clamp = True), help = 'Number of threads to use')
@click.option('-U','--unlinked', panel = "Parameters", is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples', panel = "Parameters", is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for phasing rather than those found in the inputs')
@click.help_option('--help', panel = "Workflow Options", hidden = True)
@click.argument('reference',required = True, type=FASTAfile(), nargs = 1)
@click.argument('vcf', required = True, type = VCFfile(), nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def bam(vcf, inputs, output_dir, threads, unlinked, vcf_samples, molecule_distance, reference, snakemake, extra_params, quiet, hpc, clean, container, setup):
    """
    Phase alignments using phased SNPs

    Provide the reference fasta, then a **phased** vcf file followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/myotis*.bam`), or both.
    
    Presence and type of linked-read data is auto-detected, but you may choose to omit barcode
    information with `-U`. Use `--vcf-samples` to phase only the samples present in your input
    `VCF` file rather than all the samples present in the `INPUTS` alignments.
    """
    workflow = Workflow("phase_bam", "phase_bam.smk", output_dir, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.conda = ["phase"]

    ## checks and validations ##
    alignments = XAM(inputs, detect_bc= not unlinked, quiet = quiet > 0)
    vcffile = VCF(vcf, workflow.workflow_directory, quiet = quiet > 0)
    vcffile.check_phase()
    vcffile.match_samples(alignments.files, vcf_samples)
    fasta = FASTA(reference)

    workflow.linkedreads["type"] = alignments.lr_type
    workflow.input(vcffile.file, "vcf")
    workflow.input(fasta.file, "reference")
    workflow.input(alignments.files, "alignments")
    workflow.param(molecule_distance, "distance-threshold")
    workflow.param(alignments.lr_type != "none", "use-linked-info")
    workflow.param(vcf_samples, "prioritize-vcf-samples")
    if extra_params:
        workflow.param(extra_params, "extra")

    workflow.start_text = workflow_info(
        ("Input VCF:", os.path.basename(vcffile.file)),
        ("Samples:", min(len(vcffile.samples), alignments.count)),
        ("Barcode Type:", alignments.lr_type),
        ("Reference:", os.path.basename(reference)) if reference else None,
        ("Output Folder:", os.path.relpath(output_dir) + "/")
    )

    workflow.initialize(setup)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/phase")
@click.option('-x', '--extra-params', panel = "Parameters", type = HapCutParams(), help = 'Additional HapCut2 parameters, in quotes')
@click.option('-r', '--reference', panel = "Parameters", type=FASTAfile(), help = 'Path to reference genome if wanting to also extract reads spanning indels')
@click.option('-q', '--min-map-quality', panel = "Parameters", default = 20, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality for phasing')
@click.option('-m', '--min-base-quality', panel = "Parameters", default = 13, show_default = True, type = click.IntRange(0, 100, clamp = True), help = 'Minimum base quality for phasing')
@click.option('-d', '--molecule-distance', panel = "Parameters", default = 100000, show_default = True, type = click.IntRange(min = 100), help = 'Distance cutoff to split molecules (bp)')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Phase/snp", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prune-threshold', panel = "Parameters", default = 30, show_default = True, type = click.IntRange(0,100, clamp = True), help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(2, 999, clamp = True), help = 'Number of threads to use')
@click.option('-U','--unlinked', panel = "Parameters", is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--contigs', panel = "Workflow Options",  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples', panel = "Parameters", is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for phasing rather than those found in the inputs')
@click.help_option('--help', panel = "Workflow Options", hidden = True)
@click.argument('vcf', required = True, type = VCFfile(), nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def snp(vcf, inputs, output_dir, threads, unlinked, min_map_quality, min_base_quality, molecule_distance, prune_threshold, vcf_samples, reference, snakemake, extra_params, skip_reports, quiet, hpc, clean, container, contigs, setup):
    """
    Phase SNPs into haplotypes

    Provide the vcf file followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/myotis*.bam`), or both.
    
    Presence and type of linked-read data is auto-detected, but you may choose to omit barcode
    information with `-U`. Use `--vcf-samples` to phase only the samples present in your input
    `VCF` file rather than all the samples present in the `INPUTS` alignments.
    """
    workflow = Workflow("phase_snp", "phase_snp.smk", output_dir, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.notebook_files = ["hapcut.ipynb"]
    workflow.conda = ["phase", "report"]

    ## checks and validations ##
    alignments = XAM(inputs, detect_bc= not unlinked, quiet = quiet > 0)
    vcffile = VCF(vcf, workflow.workflow_directory, quiet = quiet > 0)
    vcffile.match_samples(alignments.files, vcf_samples)
    if contigs:
        vcffile.match_contigs(contigs)

    workflow.linkedreads["type"] = alignments.lr_type
    workflow.notebooks["skip"] = skip_reports
    workflow.notebooks["plot-contigs"] = contigs if contigs else "default"
    workflow.input(vcffile.file, "vcf")
    if reference:
        fasta = FASTA(reference)
        workflow.input(fasta.file, "reference")
    workflow.input(alignments.files, "alignments")
    workflow.param(prune_threshold, "prune")
    workflow.param(min_map_quality, "min-map-quality")
    workflow.param(min_base_quality, "min-base-quality")
    workflow.param(molecule_distance, "distance-threshold")
    workflow.param(vcf_samples, "prioritize-vcf-samples")
    if extra_params:
        workflow.param(extra_params, "extra")

    workflow.start_text = workflow_info(
        ("Input VCF:", os.path.basename(vcffile.file)),
        ("Samples:", min(len(vcffile.samples), alignments.count)),
        ("Barcode Type:", alignments.lr_type),
        ("Phase Indels:", "yes" if reference else "no"),
        ("Reference:", os.path.basename(reference)) if reference else None,
        ("Output Folder:", os.path.relpath(output_dir) + "/")
    )

    workflow.initialize(setup)

phase.add_command(bam)
phase.add_command(snp)
