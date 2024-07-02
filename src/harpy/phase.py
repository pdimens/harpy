"""Harpy haplotype phasing workflow"""

import os
import sys
import subprocess
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, fetch_report
from .fileparsers import parse_alignment_inputs
from .printfunctions import print_onstart
from .validations import vcf_samplematch, validate_bam_RG, validate_input_by_ext

docstring = {
        "harpy phase": [
        {
            "name": "Parameters",
            "options": ["--vcf", "--molecule-distance", "--genome", "--prune-threshold", "--ignore-bx", "--extra-params", "--vcf-samples"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },     
    ]
}

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/phase")
@click.option('-v', '--vcf', required = True, type=click.Path(exists=True, dir_okay=False), help = 'Path to BCF/VCF file')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Base-pair distance threshold to separate molecules')
@click.option('-p', '--prune-threshold', default = 7, show_default = True, type = click.IntRange(0,100), help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-b', '--ignore-bx',  is_flag = True, show_default = True, default = False, help = 'Ignore barcodes when phasing')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for phasing rather than those found the inputs')
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), help = 'Path to genome assembly if wanting to also extract reads spanning indels')
@click.option('-x', '--extra-params', type = str, help = 'Additional HapCut2 parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 2, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Phase", show_default=True, help = 'Output directory name')
@click.option('--snakemake',  type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--config-only',  is_flag = True, hidden = True, default = False, help = 'Create the config.yaml file and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def phase(inputs, output_dir, vcf, threads, molecule_distance, prune_threshold, vcf_samples, genome, snakemake, extra_params, ignore_bx, skipreports, quiet, hpc, conda, config_only):
    """
    Phase SNPs into haplotypes

    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/myotis*.bam`), or both.
    
    You may choose to omit barcode information with `--ignore-bx`, although it's usually
    better to include that information. Use `--vcf-samples` to phase only
    the samples present in your input `--vcf` file rather than all the samples present in
    the `INPUT` alignments.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/phase.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    os.makedirs(f"{workflowdir}/input", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    samplenames = vcf_samplematch(vcf, bamlist, vcf_samples)
    validate_input_by_ext(vcf, "--vcf", ["vcf", "bcf", "vcf.gz"])
    validate_bam_RG(bamlist)
    if genome:
        validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    fetch_rule(workflowdir, "phase.smk")
    fetch_report(workflowdir, "HapCut2.Rmd")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: phase\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"noBX: {ignore_bx}\n")
        config.write(f"prune: {prune_threshold/100}\n")
        config.write(f"molecule_distance: {molecule_distance}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  variantfile: {vcf}\n")
        if genome is not None:
            config.write(f"  genome: {genome}\n")
            if not os.path.exists(f"{genome}.fai"):
                subprocess.run(f"samtools faidx --fai-idx {genome}.fai {genome}".split())
        config.write("  alignments:\n")
        for i in bamlist:
            config.write(f"    - {i}\n")
    if config_only:
        sys.exit(0)

    print_onstart(
        f"Input VCF: {vcf}\nSamples in VCF: {len(samplenames)}\nAlignments Provided: {n}\nOutput Directory: {output_dir}/",
        "phase"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)