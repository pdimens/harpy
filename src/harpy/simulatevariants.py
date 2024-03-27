from .helperfunctions import generate_conda_deps, fetch_file 
from .printfunctions import print_onstart, print_error
import rich_click as click
from pathlib import Path
import subprocess
import os
import sys

@click.command(no_args_is_help = True, epilog = "This workflow can be quite technical, please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate")
#@click.option('-p', '--parameters', type = click.Path(exists=True), help = "Simulation parameter YAML file")
#@click.option('-v', '--vcf', type=click.Path(exists=True), help = 'VCF file of known variants to simulate')
@click.option('-v', '--variant-type', type = click.Choice(["snp, indel", "inversion", "cnv", "translocation"]), show_choices=True, required=True, help = "Type of variant to simulate {`snp`,`indel`,`inversion`,`cnv`,`translocation`}")
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random variants to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum variant size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum variant size (bp)")
@click.option('-r', '--ratio', type = click.FloatRange(min=0), show_choices=True, help = "Ratio for the selected variant")
#@click.option('-a', '--indel-size-alpha', type = click.FloatRange(min=0), default = 2.0, help = "Exponent Alpha for power-law-fitted size distribution")
#@click.option('-a', '--indel-size-constant', type = click.FloatRange(min=0), default = 0.5, help = "Exponent constant for power-law-fitted size distribution")
@click.option('-z', '--cnv-max-copy', type = click.IntRange(min = 1), default=10, show_default=True, help = "Maximum number of copies")
@click.option('-l', '--cnv-ratio', type = click.FloatRange(min=0), default=1, show_default=True, help = " Relative ratio of DNA gain over DNA loss")
@click.option('-c', '--centromeres', type = click.Path(exists=True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True), help = "GFF3 file of genes to avoid when simulating (requires `--snp-coding-partition` for SNPs)")
@click.option('-p', '--snp-gene-constraints', type = click.Choice(["noncoding", "coding", "2d", "4d"]), help = "How to constrain randomly simulated SNPs {`noncoding`,`coding`,`2d`,`4d`}")
@click.option('-h', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '\% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True), help = "Text file of chromosomes to avoid")
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/variants", show_default=True, help = 'Name of output directory')
@click.argument('genome', required=True, type=click.Path(exists=True), nargs=1)
@click.argument('vcf', required=False, type=click.Path(exists=True), nargs=1)
def variants(genome, output_dir, variant_type, count, min_size, max_size, ratio, cnv_max_copy, cnv_ratio, centromeres, genes, snp_gene_constraints, heterozygosity, exclude_chr, randomseed, snakemake, quiet, print_only):
    """
    Introduce variants into a genome
 
    Create a haploid genome with variants introduced into it. Use either a VCF file to simulate known variants
    or the command line options listed below to simulate random variants of the selected type. The
    limitations of the simulator (`simuG`) are such that you may simulate only one type of variant at a time,
    so you will likely need to run this module again on the resulting genome. Setting `--heterozygosity` greater
    than `0` will also create a subsampled VCF file with which you can create another simulated genome (diploid)
    by running this module again.

    ### Notes on `--ratio`
    Where possible, command line options are consolidated for different variant types and will be interpolated properly
    within `simuG`. The helpstring below organizes which options are for which variants. Some options,
    like `--ratio` mean different things for different variant types and have special-meanining values of `9999` and `0`:
    
    | variant | ratio meaning | `9999` | `0` | default |
    |:----|:----|:---|:----| :---:|
    | snp | transitions / transversions | transit. only | transv. only | 0.5 |
    | indel | insertions / deletions | insert. only | delet. only | 1.0 |
    | cnv | tandem / dispersed | tand. only | disp. only | 1.0 |
    | `--cnv-ratio` | copy gain / loss | gain only | loss only| 1.0 | 
    """
    if variant_type == "snp":
        if (snp_gene_constraints and not genes) or (genes and not snp_gene_constraints):
            print_error("The options `--genes` and `--snp-coding-partition` must be used together for SNP variants.")
    if not vcf and count == 0:
        print_error("Please provide a value for `--count` that is greater than 0.")
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/simulate-variants.smk')
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
    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)   
    _ = Path(f"{workflowdir}/input/{os.path.basename(genome)}").symlink_to(Path(genome).absolute())
    printmsg = f"Inpute Genome: {genome}\nOutput Directory: {output_dir}/\n\Variant Type:{variant_type}\n"
    genome_link = f"{workflowdir}/input/{os.path.basename(genome)}"
    _ = Path(genome_link).symlink_to(Path(vcf).absolute())
    if vcf:
        vcf_link = f"{workflowdir}/input/{os.path.basename(vcf)}"
        _ = Path(vcf_link).symlink_to(Path(vcf).absolute())
        printmsg += f"Input VCF: {vcf}\n"
    if centromeres:
        centromeres_link = f"{workflowdir}/input/{os.path.basename(centromeres)}"
        _ = Path(centromeres_link).symlink_to(Path(centromeres).absolute())
        printmsg += f"Centromere GFF: {centromeres}\n"
    if genes:
        genes_link = f"{workflowdir}/input/{os.path.basename(genes)}"
        _ = Path(genes_link).symlink_to(Path(genes).absolute())
        printmsg += f"Genes GFF: {genes}\n"
    if exclude_chr:
        exclude_link = f"{workflowdir}/input/{os.path.basename(exclude_chr)}"
        _ = Path(exclude_link).symlink_to(Path(exclude_chr).absolute())
        printmsg += f"Excluded Chromosomes: {exclude_chr}\n"
    fetch_file("simulate-variants.smk", f"{workflowdir}/")
    fetch_file("simuG.pl", f"{workflowdir}/scripts/")
    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"input_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"variant_type: {variant_type}\n")
        config.write(f"genome: {genome}\n")
        if vcf:
            config.write(f"vcf: {vcf}\n")
        else:
            config.write(f"count: {count}\n")
            if variant_type == "snp":
                config.write(f"snp_gene_constraints: {snp_gene_constraints}\n") if snp_gene_constraints else None
                config.write(f"ratio: {ratio}\n") if ratio else None
            elif variant_type == "indel":
                config.write(f"ratio: {ratio}\n") if ratio else None
            elif variant_type == "inversion":
                config.write(f"min_size: {min_size}\n") if min_size else None
                config.write(f"max_size: {max_size}\n") if max_size else None
                config.write(f"ratio: {ratio}\n") if ratio else None
            elif variant_type == "cnv":
                config.write(f"min_size: {min_size}\n") if min_size else None
                config.write(f"max_size: {max_size}\n") if max_size else None
                config.write(f"ratio: {ratio}\n") if ratio else None
                config.write(f"cnv_max_copy: {cnv_max_copy}\n") if cnv_max_copy else None
                config.write(f"cnv_ratio: {cnv_gain_loss_ratio}\n") if cnv_ratio else None
            config.write(f"centromeres: {centromeres}\n") if centromeres else None
            config.write(f"genes: {genes}\n") if genes else None
            config.write(f"heterozygosity: {heterozygosity}\n") if heterozygosity else None
            config.write(f"exclude_chr: {exclude_chr}\n") if exclude_chr else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        printmsg.rstrip("\n"),
        "simulate variants: {variant_type}"
    )
    return 0