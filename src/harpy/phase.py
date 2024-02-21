from .helperfunctions import fetch_file, generate_conda_deps, getnames, print_onstart
from .helperfunctions import vcfcheck, vcf_samplematch, validate_bamfiles, parse_alignment_inputs
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
@click.option('-s', '--snakemake',  type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def phase(input, vcf, threads, molecule_distance, prune_threshold, vcf_samples, genome, snakemake, extra_params, ignore_bx, skipreports, quiet, print_only):
    """
    Phase SNPs into haplotypes

    You may choose to omit barcode information with `--ignore-bx`, although it's usually
    better to include that information. Use `--vcf-samples` to phase only
    the samples present in your input `--vcf` file rather than all the samples present in
    the `--directory`.
    """
    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append('Phase/workflow/phase-pop.smk')
    command.append("--configfile")
    command.append("Phase/workflow/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)

    fetch_file("phase-pop.smk", "Phase/workflow/")
    fetch_file("HapCut2.Rmd", "Phase/workflow/report/")
    vcfcheck(vcf)
    sn = parse_alignment_inputs(input, "Phase/workflow/input")
    samplenames = vcf_samplematch(vcf, "Phase/workflow/input", vcf_samples)
    validate_bamfiles("Phase/workflow/input", samplenames)
    prune_threshold /= 100

    with open("Phase/workflow/config.yml", "w") as config:
        config.write(f"seq_directory: Phase/workflow/input\n")
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
        f"Input VCF: {vcf}\nSamples in VCF: {len(samplenames)}\nOutput Directory: Phase/",
        "phase"
    )
    generate_conda_deps()
    _module = subprocess.run(command)
    sys.exit(_module.returncode)