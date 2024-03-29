from .helperfunctions import generate_conda_deps, fetch_file 
from .printfunctions import print_onstart, print_error
import rich_click as click
from pathlib import Path
import subprocess
import os
import sys

def symlink(original, destination):
    """Create a symbolic link from original -> destination if the destinationd doesn't already exist."""
    if not (Path(destination).is_symlink() or Path(destination).exists()):
        Path(destination).symlink_to(Path(original).absolute()) 

@click.command(no_args_is_help = True, epilog = "This workflow can be quite technical, please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate")
@click.option('-v', '--snp-vcf', type=click.Path(exists=True), help = 'VCF file of known snps to simulate')
@click.option('-i', '--indel-vcf', type=click.Path(exists=True), help = 'VCF file of known indels to simulate')
@click.option('-n', '--snp-count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random snps to simluate")
@click.option('-m', '--indel-count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random indels to simluate")
@click.option('-r', '--titv-ratio', type = click.FloatRange(min=0), default= 0.5, show_choices=True, show_default=True, help = "Transition/Transversion ratio for snps")
@click.option('-d', '--indel-ratio', type = click.FloatRange(min=0), default= 1, show_choices=True, show_default=True, help = "Insertion/Deletion ratio for indels")
@click.option('-a', '--indel-size-alpha', type = click.FloatRange(min=0), hidden = True, default = 2.0, help = "Exponent Alpha for power-law-fitted size distribution")
@click.option('-l', '--indel-size-constant', type = click.FloatRange(min=0), default = 0.5, hidden = True, help = "Exponent constant for power-law-fitted size distribution")
@click.option('-c', '--centromeres', type = click.Path(exists=True), help = "GFF3 file of centromeres to avoid")
@click.option('-z', '--snp-gene-constraints', type = click.Choice(["noncoding", "coding", "2d", "4d"]), help = "How to constrain randomly simulated SNPs {`noncoding`,`coding`,`2d`,`4d`}")
@click.option('-g', '--genes', type = click.Path(exists=True), help = "GFF3 file of genes to use with `--snp-gene-constraints`")
@click.option('-h', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '\% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True), help = "Text file of chromosomes to avoid")
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('-p', '--prefix', type = str, default= "simulation", show_default=True, help = "Naming prefix for output files")
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/snpindel", show_default=True, help = 'Name of output directory')
@click.argument('genome', required=True, type=click.Path(exists=True), nargs=1)
def snpindel(genome, snp_vcf, indel_vcf, output_dir, prefix, snp_count, indel_count, titv_ratio, indel_ratio, indel_size_alpha, indel_size_constant, centromeres, genes, snp_gene_constraints, heterozygosity, exclude_chr, randomseed, snakemake, quiet, print_only):
    """
    Introduce snps and/or indels into a genome
 
    Create a haploid genome with snps and/or indels introduced into it. Use either a VCF file to simulate known variants
    via `--snp-vcf` and/or `--indel-vcf` or the command line options listed below to simulate random variants of the selected
    type. Setting `--heterozygosity` greater than `0` will also create a subsampled VCF file with which you can create
    another simulated genome (diploid) by running the `simulate` module again.
    The ratio parameters control different things for snp and indel variants and have special meanings when setting
    the value to either `9999` or `0` :
    
    | ratio | meaning | `9999` | `0` |
    |:---- |:---- |:--- |:---- |
    | `--snp-ratio`   | transitions / transversions | transit. only | transv. only |
    | `--indel-ratio` | insertions / deletions      | insert. only  | delet. only  |
    """
    if (snp_gene_constraints and not genes) or (genes and not snp_gene_constraints):
        print_error("The options `--genes` and `--snp-coding-partition` must be used together for SNP variants.")
    if (not snp_vcf and snp_count == 0) or (not indel_vcf and indel_count == 0):
        print_error("You must either provide a vcf file of known variants to simulate or a count of that variant to randomly simulate.")
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores 1 --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/simulate-snpindel.smk')
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
    printmsg = f"Inpute Genome: {genome}\nOutput Directory: {output_dir}/\n"
    genome_link = f"{workflowdir}/input/{os.path.basename(genome)}"
    _ = Path(genome_link).symlink_to(Path(genome).absolute())
    if snp_vcf:
        snp_vcf_link = f"{workflowdir}/input/{os.path.basename(snp_vcf)}"
        _ = Path(snp_vcf_link).symlink_to(Path(snp_vcf).absolute())
        printmsg += f"Input vcf (snp): {vcf}\n"
    if indel_vcf:
        indel_vcf_link = f"{workflowdir}/input/{os.path.basename(indel_vcf)}"
        _ = Path(indel_vcf_link).symlink_to(Path(indel_vcf).absolute())
        printmsg += f"Input vcf (indel): {vcf}\n"
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
        config.write(f"genome: {genome}\n")
        if snp_vcf:
            config.write(f"snp_vcf: {snp_vcf}\n")
        if indel_vcf:
            config.write(f"indel_vcf: {indel_vcf}\n")
        else:
            config.write(f"snp_count: {snp_count}\n") if snp_count else None
            config.write(f"indel_count: {snp_count}\n") if indel_count else None
            config.write(f"snp_gene_constraints: {snp_gene_constraints}\n") if snp_gene_constraints else None
            config.write(f"titv_ratio: {titv_ratio}\n") if ratio else None
            config.write(f"indel_ratio: {indel_ratio}\n") if ratio else None
            config.write(f"indel_size_alpha: {indel_size_alpha}\n") if indel_size_alpha else None
            config.write(f"indel_size_constant: {indel_size_constant}\n") if indel_size_constant else None
            config.write(f"centromeres: {centromeres}\n") if centromeres else None
            config.write(f"genes: {genes}\n") if genes else None
            config.write(f"heterozygosity: {heterozygosity}\n") if heterozygosity else None
            config.write(f"exclude_chr: {exclude_chr}\n") if exclude_chr else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        printmsg.rstrip("\n"),
        "simulate variants: snpindel"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "Please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate")
@click.option('-v', '--vcf', type=click.Path(exists=True), help = 'VCF file of known inversions to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random inversions to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum inversion size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum inversion size (bp)")
@click.option('-c', '--centromeres', type = click.Path(exists=True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-h', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '\% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True), help = "Text file of chromosomes to avoid")
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('-p', '--prefix', type = str, default= "simulation", show_default=True, help = "Naming prefix for output files")
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/inversion", show_default=True, help = 'Name of output directory')
@click.argument('genome', required=True, type=click.Path(exists=True), nargs=1)
def inversion(genome, vcf, prefix, output_dir, count, min_size, max_size, centromeres, genes, heterozygosity, exclude_chr, randomseed, snakemake, quiet, print_only):
    """
    Introduce inversions into a genome
 
    Create a haploid genome with inversions introduced into it. Use either a VCF file to simulate known inversions
    or the command line options listed below to simulate random inversions. Setting `--heterozygosity` greater
    than `0` will also create a subsampled VCF file with which you can create another simulated genome (diploid)
    by running this module again.
    """
    if not vcf and count == 0:
        print_error("Provide either a `--count` of cnv to randomly simulate or a `--vcf` of known variants to simulate.")

    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores 1 --directory .'.split()
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
    printmsg = f"Inpute Genome: {genome}\nOutput Directory: {output_dir}/\n"
    genome_link = f"{workflowdir}/input/{os.path.basename(genome)}"
    _ = Path(genome_link).symlink_to(Path(genome).absolute())
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
            config.write(f"min_size: {min_size}\n") if min_size else None
            config.write(f"max_size: {max_size}\n") if max_size else None
            config.write(f"centromeres: {centromeres}\n") if centromeres else None
            config.write(f"genes: {genes}\n") if genes else None
            config.write(f"heterozygosity: {heterozygosity}\n") if heterozygosity else None
            config.write(f"exclude_chr: {exclude_chr}\n") if exclude_chr else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        printmsg.rstrip("\n"),
        "simulate variants: inversion"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "Please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate")
@click.option('-v', '--vcf', type=click.Path(exists=True), help = 'VCF file of known copy number variants to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random variants to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum variant size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum variant size (bp)")
@click.option('-d', '--dup-ratio', type = click.FloatRange(min=0), default = 1, show_choices=True, show_default=True, help = "Ratio for the selected variant")
@click.option('-l', '--gain-ratio', type = click.FloatRange(min=0), default=1, show_default=True, help = " Relative ratio of DNA gain over DNA loss")
@click.option('-z', '--max-copy', type = click.IntRange(min = 1), default=10, show_default=True, help = "Maximum number of copies")
@click.option('-c', '--centromeres', type = click.Path(exists=True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True), help = "GFF3 file of genes to avoid when simulating (requires `--snp-coding-partition` for SNPs)")
@click.option('-h', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '\% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True), help = "Text file of chromosomes to avoid")
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('-p', '--prefix', type = str, default= "simulate.cnv", show_default=True, help = "Naming prefix for output files")
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/cnv", show_default=True, help = 'Name of output directory')
@click.argument('genome', required=True, type=click.Path(exists=True), nargs=1)
def cnv(genome, output_dir, vcf, prefix, count, min_size, max_size, dup_ratio, max_copy, gain_ratio, centromeres, genes, heterozygosity, exclude_chr, randomseed, snakemake, quiet, print_only):
    """
    Introduce copy number variants into a genome
 
    Create a haploid genome with copy number variants introduced into it. Use either a VCF file to simulate known CNVs
    or the command line options listed below to simulate random variants of the selected type. Setting `--heterozygosity` greater
    than `0` will also create a subsampled VCF file with which you can create another simulated genome (diploid)
    by running this module again.

    The two ratio parameters control different things and have special meanings when setting their value to either `9999` or `0`:
    
    | ratio | meaning | `9999` | `0` |
    |:----|:----|:---|:----|
    | `--dup-ratio` | tandem / dispersed | tand. only | disp. only |
    | `--gain-ratio` | copy gain / loss | gain only | loss only| 
    """
    if not vcf and count == 0:
        print_error("Provide either a `--count` of cnv to randomly simulate or a `--vcf` of known cnv to simulate.")
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores 1 --directory .'.split()
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
    symlink(genome, f"{workflowdir}/input/{os.path.basename(genome)}")
    printmsg = f"Inpute Genome: {genome}\nOutput Directory: {output_dir}/\n"
    genome_link = f"{workflowdir}/input/{os.path.basename(genome)}"
    symlink(genome, genome_link)
    if vcf:
        vcf_link = f"{workflowdir}/input/{os.path.basename(vcf)}"
        symlink(vcf, vcf_link)
        printmsg += f"Input VCF: {vcf}\n"
    else:
        printmsg += f"Mode: Random variants\n"
    if centromeres:
        centromeres_link = f"{workflowdir}/input/{os.path.basename(centromeres)}"
        symlink(centromeres, centromeres_link)
        printmsg += f"Centromere GFF: {centromeres}\n"
    if genes:
        genes_link = f"{workflowdir}/input/{os.path.basename(genes)}"
        symlink(genes, genes_link)
        printmsg += f"Genes GFF: {genes}\n"
    if exclude_chr:
        exclude_link = f"{workflowdir}/input/{os.path.basename(exclude_chr)}"
        symlink(exclude_chr, exclude_link)
        printmsg += f"Excluded Chromosomes: {exclude_chr}\n"
    fetch_file("simulate-variants.smk", f"{workflowdir}/")
    fetch_file("simuG.pl", f"{workflowdir}/scripts/")
    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"input_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"variant_type: cnv\n")
        config.write(f"genome: {genome_link}\n")
        config.write(f"prefix: {prefix}\n")
        if vcf:
            config.write(f"vcf: {vcf_link}\n")
        else:
            config.write(f"count: {count}\n")
            config.write(f"min_size: {min_size}\n") if min_size else None
            config.write(f"max_size: {max_size}\n") if max_size else None
            config.write(f"dup_ratio: {dup_ratio}\n") if dup_ratio else None
            config.write(f"cnv_max_copy: {max_copy}\n") if max_copy else None
            config.write(f"gain_ratio: {gain_ratio}\n") if gain_ratio else None
            config.write(f"centromeres: {centromeres_link}\n") if centromeres else None
            config.write(f"genes: {genes_link}\n") if genes else None
            config.write(f"exclude_chr: {exclude_link}\n") if exclude_chr else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"heterozygosity: {heterozygosity}\n")
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        printmsg.rstrip("\n"),
        "simulate cnv"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "Please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate")
@click.option('-v', '--vcf', type=click.Path(exists=True), help = 'VCF file of known translocations to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random translocations to simluate")
@click.option('-c', '--centromeres', type = click.Path(exists=True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-h', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '\% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True), help = "Text file of chromosomes to avoid")
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('-p', '--prefix', type = str, default= "simulation", show_default=True, help = "Naming prefix for output files")
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/translocation", show_default=True, help = 'Name of output directory')
@click.argument('genome', required=True, type=click.Path(exists=True), nargs=1)
def translocation(genome, output_dir, prefix, vcf, count, centromeres, genes, snp_gene_constraints, heterozygosity, exclude_chr, randomseed, snakemake, quiet, print_only):
    """
    Introduce transolcations into a genome
 
    Create a haploid genome with translocations introduced into it. Use either a VCF file to 
    simulate known translocations or the command line options listed below to simulate random
    variants of the selected type. Setting `--heterozygosity` greater than `0` will also create
    a subsampled VCF file with which you can create another simulated genome (diploid) by running
    this module again.
    """
    if not vcf and count == 0:
        print_error("Provide either a `--count` of cnv to randomly simulate or a `--vcf` of known cnv to simulate.")
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores 1 --directory .'.split()
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
    printmsg = f"Inpute Genome: {genome}\nOutput Directory: {output_dir}/\n"
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
            config.write(f"centromeres: {centromeres}\n") if centromeres else None
            config.write(f"genes: {genes}\n") if genes else None
            config.write(f"heterozygosity: {heterozygosity}\n") if heterozygosity else None
            config.write(f"exclude_chr: {exclude_chr}\n") if exclude_chr else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        printmsg.rstrip("\n"),
        "simulate variants: translocation"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)


#TODO ONLY CNV WAS FIXED. OTHERS NEED FIXING