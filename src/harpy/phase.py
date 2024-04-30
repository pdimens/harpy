from .helperfunctions import fetch_rule, fetch_report, generate_conda_deps
from .fileparsers import parse_alignment_inputs
from .printfunctions import print_onstart
from .validations import vcfcheck, vcf_samplematch, validate_bamfiles, validate_input_by_ext
import sys
import os
import subprocess
import rich_click as click

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/phase")
@click.option('-v', '--vcf', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Path to BCF/VCF file')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance threshold to separate molecules')
@click.option('-p', '--prune-threshold', default = 7, show_default = True, type = click.IntRange(0,100), metavar = "Integer", help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-b', '--ignore-bx',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Ignore barcodes when phasing')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Use samples present in vcf file for phasing rather than those found the directory')
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), metavar = "File path", help = 'Path to genome assembly if wanting to also extract reads spanning indels')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional HapCut2 parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 2, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Phase", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('--snakemake',  type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def phase(input, output_dir, vcf, threads, molecule_distance, prune_threshold, vcf_samples, genome, snakemake, extra_params, ignore_bx, skipreports, quiet, print_only):
    """
    Phase SNPs into haplotypes

    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/myotis*.bam`), or both.
    
    You may choose to omit barcode information with `--ignore-bx`, although it's usually
    better to include that information. Use `--vcf-samples` to phase only
    the samples present in your input `--vcf` file rather than all the samples present in
    the input alignments.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f"{workflowdir}/phase.smk")
    command.append("--configfile")
    command.append(f"{workflowdir}/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit(0)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    sn = parse_alignment_inputs(input, f"{workflowdir}/input")
    samplenames = vcf_samplematch(vcf, f"{workflowdir}/input", vcf_samples)
    vcfcheck(vcf)
    validate_bamfiles(f"{workflowdir}/input", samplenames)
    if genome:
        validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    fetch_rule(workflowdir, "phase.smk")
    fetch_report(workflowdir, "HapCut2.Rmd")
    prune_threshold /= 100

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplenames: {samplenames}\n")
        config.write(f"variantfile: {vcf}\n")
        config.write(f"noBX: {ignore_bx}\n")
        config.write(f"prune: {prune_threshold}\n")
        config.write(f"molecule_distance: {molecule_distance}\n")
        if genome is not None:
            config.write(f"indels: {genome}\n")
            if not os.path.exists(f"{genome}.fai"):
                subprocess.run(f"samtools faidx --fai-idx {genome}.fai {genome}".split())
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Input VCF: {vcf}\nSamples in VCF: {len(samplenames)}\nAlignments Provided: {len(sn)}\nOutput Directory: {output_dir}/",
        "phase"
    )
    return command
