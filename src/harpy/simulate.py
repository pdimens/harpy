"""Harpy workflows to simulate genomic variants and linked-reads"""
#TODO REPLACE SYMLINK CODE WITH HPC CONDITIONAL
import os
import sys
import shutil
import subprocess
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, fetch_script, symlink
from .printfunctions import print_onstart, print_error
from .validations import validate_input_by_ext

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def simulate():
    """
    Simulate variants or linked reads from a genome

    To simulate genomic variants, provide an additional subcommand {`snpindel`,`inversion`,`cnv`,`translocation`} 
    to get more information about that workflow. The limitations of the simulator
    (`simuG`) are such that you may simulate only one type of variant at a time,
    so you may need to run this module again on the resulting genome. Use `simulate linkedreads`
    to simulate haplotag linked-reads from a diploid genome, which you can create by simulating
    genomic variants.
    """

commandstring = {
    "harpy simulate": [
        {
            "name": "Linked Read Sequences",
            "commands": ["linkedreads"],
        },
        {
            "name": "Genomic Variants",
            "commands": ["snpindel","inversion", "cnv", "translocation"],
        }
    ]
}

docstring = {
    "harpy simulate linkedreads": [
        {
            "name": "Parameters",
            "options": ["--barcodes", "--read-pairs", "--outer-distance", "--distance-sd", "--mutation-rate", "--molecule-length", "--partitions", "--molecules-per"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },     
    ],
    "harpy simulate snpindel": [
        {
            "name": "Known Variants",
            "options": ["--snp-vcf", "--indel-vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--snp-count", "--indel-count", "--titv-ratio", "--indel-ratio", "--snp-gene-constraints", "--genes", "--centromeres", "--exclude-chr"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--prefix", "--heterozygosity", "--randomseed", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy simulate inversion": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--count", "--min-size", "--max-size", "--genes", "--centromeres", "--exclude-chr"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--prefix", "--heterozygosity", "--randomseed", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy simulate cnv": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--count", "--min-size", "--max-size", "--max-copy", "--dup-ratio", "--gain-ratio", "--genes", "--centromeres", "--exclude-chr"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--prefix", "--heterozygosity", "--randomseed", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy simulate translocation": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--count", "--genes", "--centromeres", "--exclude-chr"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir","--prefix","--heterozygosity", "--randomseed", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-linkedreads")
@click.option('-d', '--outer-distance', type = click.IntRange(min = 100), default = 350, show_default= True, help = "Outer distance between paired-end reads (bp)")
@click.option('-s', '--distance-sd', type = click.IntRange(min = 1), default = 15, show_default=True,  help = "Standard deviation of read-pair distance")
@click.option('-b', '--barcodes', type = click.Path(exists=True, dir_okay=False), help = "File of linked-read barcodes to add to reads")
@click.option('-n', '--read-pairs', type = click.FloatRange(min = 0.001), default = 600, show_default=True,  help = "Number (in millions) of read pairs to simulate")
@click.option('-l', '--molecule-length', type = click.IntRange(min = 2), default = 100, show_default=True,  help = "Mean molecule length (kbp)")
@click.option('-r', '--mutation-rate', type = click.FloatRange(min = 0), default=0.001, show_default=True,  help = "Random mutation rate for simulating reads")
@click.option('-p', '--partitions', type = click.IntRange(min = 1), default=1500, show_default=True,  help = "Number (in thousands) of partitions/beads to generate")
@click.option('-m', '--molecules-per', type = click.IntRange(min = 1), default = 10, show_default=True,  help = "Average number of molecules per partition")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/linkedreads", help = 'Name of output directory')
@click.argument('genome_hap1', required=True, type=click.Path(exists=True, dir_okay=False), nargs=1)
@click.argument('genome_hap2', required=True, type=click.Path(exists=True, dir_okay=False), nargs=1)
def linkedreads(genome_hap1, genome_hap2, output_dir, outer_distance, mutation_rate, distance_sd, barcodes, read_pairs, molecule_length, partitions, molecules_per, threads, snakemake, quiet, hpc, conda, print_only):
    """
    Create linked reads from a genome
 
    Two haplotype genomes (un/compressed fasta) need to be provided as inputs at the end of the command. If
    you don't have a diploid genome, you can simulate one with `harpy simulate` as described [in the documentation](https://pdimens.github.io/harpy/modules/simulate/simulate-variants/#simulate-diploid-assembly).

    If not providing a text file of `--barcodes`, Harpy will download the `4M-with-alts-february-2016.txt`
    file containing the standard 16-basepair 10X barcodes, which is available from 10X genomics and the
    LRSIM [GitHub repository](https://github.com/aquaskyline/LRSIM/). Barcodes in the `--barcodes` file
    are expected to be one 16-basepar barcode per line.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock  --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/simulate-linkedreads.smk "
    command += f"--configfile {workflowdir}/config.yml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake
    if print_only:
        click.echo(command)
        sys.exit(0)

    validate_input_by_ext(genome_hap1, "GENOME_HAP1", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    validate_input_by_ext(genome_hap2, "GENOME_HAP2", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fetch_rule(workflowdir, "simulate-linkedreads.smk")
    fetch_script(workflowdir, "LRSIMharpy.pl")

    with open(f"{workflowdir}/config.yml", "w", encoding="utf-8") as config:
        config.write(f"genome_hap1: {genome_hap1}\n")
        config.write(f"genome_hap2: {genome_hap2}\n")
        config.write(f"output_directory: {output_dir}\n")
        if barcodes:
            config.write(f"barcodes: {barcodes}\n")
        config.write(f"outer_distance: {outer_distance}\n")
        config.write(f"distance_sd: {distance_sd}\n")
        config.write(f"read_pairs: {read_pairs}\n")
        config.write(f"mutation_rate: {mutation_rate}\n")
        config.write(f"molecule_length: {molecule_length}\n")
        config.write(f"partitions: {partitions}\n")
        config.write(f"molecules_per_partition: {molecules_per}\n")
        config.write(f"workflow_call: {command}\n")

    onstart_text = f"Genome Haplotype 1: {os.path.basename(genome_hap1)}\n"
    onstart_text += f"Genome Haplotype 2: {os.path.basename(genome_hap2)}\n"
    onstart_text += f"Barcodes: {os.path.basename(barcodes)}\n" if barcodes else "Barcodes: 10X Default\n"
    onstart_text += f"Output Directory: {output_dir}/"
    print_onstart(onstart_text, "simulate reads")
    
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "This workflow can be quite technical, please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-variants")
@click.option('-s', '--snp-vcf', type=click.Path(exists=True, dir_okay=False), help = 'VCF file of known snps to simulate')
@click.option('-i', '--indel-vcf', type=click.Path(exists=True, dir_okay=False), help = 'VCF file of known indels to simulate')
@click.option('-n', '--snp-count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random snps to simluate")
@click.option('-m', '--indel-count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random indels to simluate")
@click.option('-r', '--titv-ratio', type = click.FloatRange(min=0), default= 0.5, show_choices=True, show_default=True, help = "Transition/Transversion ratio for snps")
@click.option('-d', '--indel-ratio', type = click.FloatRange(min=0), default= 1, show_choices=True, show_default=True, help = "Insertion/Deletion ratio for indels")
@click.option('-a', '--indel-size-alpha', type = click.FloatRange(min=0), hidden = True, default = 2.0, help = "Exponent Alpha for power-law-fitted size distribution")
@click.option('-l', '--indel-size-constant', type = click.FloatRange(min=0), default = 0.5, hidden = True, help = "Exponent constant for power-law-fitted size distribution")
@click.option('-c', '--centromeres', type = click.Path(exists=True, dir_okay=False), help = "GFF3 file of centromeres to avoid")
@click.option('-y', '--snp-gene-constraints', type = click.Choice(["noncoding", "coding", "2d", "4d"]), help = "How to constrain randomly simulated SNPs {`noncoding`,`coding`,`2d`,`4d`}")
@click.option('-g', '--genes', type = click.Path(exists=True), help = "GFF3 file of genes to use with `--snp-gene-constraints`")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = str, default = "Simulate/snpindel", show_default=True, help = 'Name of output directory')
@click.option('-p', '--prefix', type = str, default= "sim.snpindel", show_default=True, help = "Naming prefix for output files")
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('genome', required=True, type=click.Path(exists=True, dir_okay=False), nargs=1)
def snpindel(genome, snp_vcf, indel_vcf, output_dir, prefix, snp_count, indel_count, titv_ratio, indel_ratio, indel_size_alpha, indel_size_constant, centromeres, genes, snp_gene_constraints, heterozygosity, exclude_chr, randomseed, snakemake, quiet, hpc, conda, print_only):
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
    if (not snp_vcf and snp_count == 0) and (not indel_vcf and indel_count == 0):
        print_error("You must either provide a vcf file of known variants to simulate or a count of that variant to randomly simulate.")
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores 1 --directory . '
    command += f"--snakefile {workflowdir}/simulate-snpindel.smk "
    command += f"--configfile {workflowdir}/config.yml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake
    if print_only:
        click.echo(command)
        sys.exit(0)
    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)   
    validate_input_by_ext(genome, "GENOME", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    genome_link = f"{workflowdir}/input/{os.path.basename(genome)}"
    if hpc:
        shutil.copy(genome, genome_link)
    else:
        symlink(genome, genome_link)
    printmsg = f"Inpute Genome: {genome}\nOutput Directory: {output_dir}/\n"
    if snp_vcf:
        validate_input_by_ext(snp_vcf, "--snp-vcf", ["vcf","vcf.gz","bcf"])
        snp_vcf_link = f"{workflowdir}/input/{os.path.basename(snp_vcf)}"
        symlink(snp_vcf, snp_vcf_link)
        printmsg += f"SNPs: from vcf ({snp_vcf})\n"
    elif snp_count > 0:
        printmsg += "SNPs: random\n"
    if indel_vcf:
        validate_input_by_ext(indel_vcf, "--indel-vcf", ["vcf","vcf.gz","bcf"])
        indel_vcf_link = f"{workflowdir}/input/{os.path.basename(indel_vcf)}"
        symlink(indel_vcf, indel_vcf_link)
        printmsg += f"Indels: from vcf: ({indel_vcf})\n"
    elif indel_count > 0:
        printmsg += "Indels: random\n"
    if centromeres:
        validate_input_by_ext(centromeres, "--centromeres", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        centromeres_link = f"{workflowdir}/input/{os.path.basename(centromeres)}"
        symlink(centromeres, centromeres_link)
        printmsg += f"Centromere GFF: {centromeres}\n"
    if genes:
        validate_input_by_ext(genes, "--genes", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        genes_link = f"{workflowdir}/input/{os.path.basename(genes)}"
        symlink(genes, genes_link)
        printmsg += f"Genes GFF: {genes}\n"
    if exclude_chr:
        exclude_link = f"{workflowdir}/input/{os.path.basename(exclude_chr)}"
        symlink(exclude_chr, exclude_link)
        printmsg += f"Excluded Chromosomes: {exclude_chr}\n"
    fetch_rule(workflowdir, "simulate-snpindel.smk")
    fetch_script(workflowdir, "simuG.pl")
    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yml", "w", encoding="utf-8") as config:
        config.write(f"input_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"genome: {genome_link}\n")
        config.write(f"prefix: {prefix}\n")
        if snp_vcf:
            config.write(f"snp_vcf: {snp_vcf_link}\n")
        if indel_vcf:
            config.write(f"indel_vcf: {indel_vcf_link}\n")
        else:
            config.write(f"snp_count: {snp_count}\n") if snp_count else None
            config.write(f"indel_count: {snp_count}\n") if indel_count else None
            config.write(f"snp_gene_constraints: {snp_gene_constraints}\n") if snp_gene_constraints else None
            config.write(f"titv_ratio: {titv_ratio}\n") if titv_ratio else None
            config.write(f"indel_ratio: {indel_ratio}\n") if indel_ratio else None
            config.write(f"indel_size_alpha: {indel_size_alpha}\n") if indel_size_alpha else None
            config.write(f"indel_size_constant: {indel_size_constant}\n") if indel_size_constant else None
            config.write(f"centromeres: {centromeres_link}\n") if centromeres else None
            config.write(f"genes: {genes_link}\n") if genes else None
            config.write(f"exclude_chr: {exclude_link}\n") if exclude_chr else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"heterozygosity: {heterozygosity}\n")
        config.write(f"workflow_call: {command}\n")

    print_onstart(printmsg.rstrip("\n"), "simulate variants: snpindel")
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "Please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False), help = 'VCF file of known inversions to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random inversions to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum inversion size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum inversion size (bp)")
@click.option('-c', '--centromeres', type = click.Path(exists=True, dir_okay=False), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True, dir_okay=False), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False), help = "Text file of chromosomes to avoid")
@click.option('-p', '--prefix', type = str, default= "sim.inversion", show_default=True, help = "Naming prefix for output files")
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Simulate/inversion", show_default=True, help = 'Name of output directory')
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('genome', required=True, type=click.Path(exists=True, dir_okay=False), nargs=1)
def inversion(genome, vcf, prefix, output_dir, count, min_size, max_size, centromeres, genes, heterozygosity, exclude_chr, randomseed, snakemake, quiet, hpc, conda, print_only):
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
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores 1 --directory . '
    command += f"--snakefile {workflowdir}/simulate-variants.smk "
    command += f"--configfile {workflowdir}/config.yml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake
    if print_only:
        click.echo(command)
        sys.exit(0)
    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    validate_input_by_ext(genome, "GENOME", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    genome_link = f"{workflowdir}/input/{os.path.basename(genome)}"
    if hpc:
        shutil.copy(genome, genome_link)
    else:
        symlink(genome, genome_link)
    printmsg = f"Inpute Genome: {genome}\nOutput Directory: {output_dir}/\n"

    if vcf:
        validate_input_by_ext(vcf, "--vcf", ["vcf","vcf.gz","bcf"])
        vcf_link = f"{workflowdir}/input/{os.path.basename(vcf)}"
        symlink(vcf, vcf_link)
        printmsg += f"Input VCF: {vcf}\n"
    else:
        printmsg += "Mode: Random variants\n"
    if centromeres:
        validate_input_by_ext(centromeres, "--centromeres", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        centromeres_link = f"{workflowdir}/input/{os.path.basename(centromeres)}"
        symlink(centromeres, centromeres_link)
        printmsg += f"Centromere GFF: {centromeres}\n"
    if genes:
        validate_input_by_ext(genes, "--genes", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        genes_link = f"{workflowdir}/input/{os.path.basename(genes)}"
        symlink(genes, genes_link)
        printmsg += f"Genes GFF: {genes}\n"
    if exclude_chr:
        exclude_link = f"{workflowdir}/input/{os.path.basename(exclude_chr)}"
        symlink(exclude_chr, exclude_link)
        printmsg += f"Excluded Chromosomes: {exclude_chr}\n"
    fetch_rule(workflowdir, "simulate-variants.smk")
    fetch_script(workflowdir, "simuG.pl")
    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yml", "w", encoding="utf-8") as config:
        config.write(f"input_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write("variant_type: inversion\n")
        config.write(f"genome: {genome_link}\n")
        config.write(f"prefix: {prefix}\n")
        if vcf:
            config.write(f"vcf: {vcf_link}\n")
        else:
            config.write(f"count: {count}\n")
            config.write(f"min_size: {min_size}\n") if min_size else None
            config.write(f"max_size: {max_size}\n") if max_size else None
            config.write(f"centromeres: {centromeres_link}\n") if centromeres else None
            config.write(f"genes: {genes_link}\n") if genes else None
            config.write(f"exclude_chr: {exclude_link}\n") if exclude_chr else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"heterozygosity: {heterozygosity}\n")
        config.write(f"workflow_call: {command}\n")

    print_onstart(printmsg.rstrip("\n"), "simulate variants: inversion")
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "Please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False), help = 'VCF file of known copy number variants to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random variants to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum variant size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum variant size (bp)")
@click.option('-d', '--dup-ratio', type = click.FloatRange(min=0), default = 1, show_choices=True, show_default=True, help = "Ratio for the selected variant")
@click.option('-l', '--gain-ratio', type = click.FloatRange(min=0), default=1, show_default=True, help = " Relative ratio of DNA gain over DNA loss")
@click.option('-y', '--max-copy', type = click.IntRange(min = 1), default=10, show_default=True, help = "Maximum number of copies")
@click.option('-c', '--centromeres', type = click.Path(exists=True, dir_okay=False), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True, dir_okay=False), help = "GFF3 file of genes to avoid when simulating (requires `--snp-coding-partition` for SNPs)")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = str, default = "Simulate/cnv", show_default=True, help = 'Name of output directory')
@click.option('-p', '--prefix', type = str, default= "sim.cnv", show_default=True, help = "Naming prefix for output files")
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('genome', required=True, type=click.Path(exists=True, dir_okay=False), nargs=1)
def cnv(genome, output_dir, vcf, prefix, count, min_size, max_size, dup_ratio, max_copy, gain_ratio, centromeres, genes, heterozygosity, exclude_chr, randomseed, snakemake, quiet, hpc, conda, print_only):
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
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores 1 --directory . '
    command += f"--snakefile {workflowdir}/simulate-variants.smk "
    command += f"--configfile {workflowdir}/config.yml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake
    if print_only:
        click.echo(command)
        sys.exit(0)
    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    validate_input_by_ext(genome, "GENOME", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    genome_link = f"{workflowdir}/input/{os.path.basename(genome)}"
    if hpc:
        shutil.copy(genome, genome_link)
    else:
        symlink(genome, genome_link)
    printmsg = f"Inpute Genome: {genome}\nOutput Directory: {output_dir}/\n"

    if vcf:
        validate_input_by_ext(vcf, "--vcf", ["vcf","vcf.gz","bcf"])
        vcf_link = f"{workflowdir}/input/{os.path.basename(vcf)}"
        symlink(vcf, vcf_link)
        printmsg += f"Input VCF: {vcf}\n"
    else:
        printmsg += "Mode: Random variants\n"
    if centromeres:
        validate_input_by_ext(centromeres, "--centromeres", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        centromeres_link = f"{workflowdir}/input/{os.path.basename(centromeres)}"
        symlink(centromeres, centromeres_link)
        printmsg += f"Centromere GFF: {centromeres}\n"
    if genes:
        validate_input_by_ext(genes, "--genes", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        genes_link = f"{workflowdir}/input/{os.path.basename(genes)}"
        symlink(genes, genes_link)
        printmsg += f"Genes GFF: {genes}\n"
    if exclude_chr:
        exclude_link = f"{workflowdir}/input/{os.path.basename(exclude_chr)}"
        symlink(exclude_chr, exclude_link)
        printmsg += f"Excluded Chromosomes: {exclude_chr}\n"
    fetch_rule(workflowdir, "simulate-variants.smk")
    fetch_script(workflowdir, "simuG.pl")
    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yml", "w", encoding="utf-8") as config:
        config.write(f"input_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write("variant_type: cnv\n")
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
        config.write(f"workflow_call: {command}\n")

    print_onstart(printmsg.rstrip("\n"),"simulate cnv")
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "Please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False), help = 'VCF file of known translocations to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random translocations to simluate")
@click.option('-c', '--centromeres', type = click.Path(exists=True, dir_okay=False), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True, dir_okay=False), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = str, default = "Simulate/translocation", show_default=True, help = 'Name of output directory')
@click.option('-p', '--prefix', type = str, default= "sim.translocation", show_default=True, help = "Naming prefix for output files")
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('genome', required=True, type=click.Path(exists=True, dir_okay=False), nargs=1)
def translocation(genome, output_dir, prefix, vcf, count, centromeres, genes, heterozygosity, exclude_chr, randomseed, snakemake, quiet, hpc, conda, print_only):
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
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores 1 --directory . '
    command += f"--snakefile {workflowdir}/simulate-variants.smk "
    command += f"--configfile {workflowdir}/config.yml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake
    if print_only:
        click.echo(command)
        sys.exit(0)
    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    validate_input_by_ext(genome, "GENOME", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    genome_link = f"{workflowdir}/input/{os.path.basename(genome)}"
    if hpc:
        shutil.copy(genome, genome_link)
    else:
        symlink(genome, genome_link)
    printmsg = f"Inpute Genome: {genome}\nOutput Directory: {output_dir}/\n"
    
    if vcf:
        validate_input_by_ext(vcf, "--vcf", ["vcf","vcf.gz","bcf"])
        vcf_link = f"{workflowdir}/input/{os.path.basename(vcf)}"
        symlink(vcf, vcf_link)
        printmsg += f"Input VCF: {vcf}\n"
    else:
        printmsg += "Mode: Random variants\n"
    if centromeres:
        validate_input_by_ext(centromeres, "--centromeres", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        centromeres_link = f"{workflowdir}/input/{os.path.basename(centromeres)}"
        symlink(centromeres, centromeres_link)
        printmsg += f"Centromere GFF: {centromeres}\n"
    if genes:
        validate_input_by_ext(genes, "--genes", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        genes_link = f"{workflowdir}/input/{os.path.basename(genes)}"
        symlink(genes, genes_link)
        printmsg += f"Genes GFF: {genes}\n"
    if exclude_chr:
        exclude_link = f"{workflowdir}/input/{os.path.basename(exclude_chr)}"
        symlink(exclude_chr, exclude_link)
        printmsg += f"Excluded Chromosomes: {exclude_chr}\n"
    fetch_rule(workflowdir, "simulate-variants.smk")
    fetch_script(workflowdir, "simuG.pl")
    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yml", "w", encoding="utf-8") as config:
        config.write(f"input_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write("variant_type: translocation\n")
        config.write(f"genome: {genome_link}\n")
        config.write(f"prefix: {prefix}\n")
        if vcf:
            config.write(f"vcf: {vcf_link}\n")
        else:
            config.write(f"count: {count}\n")
            config.write(f"centromeres: {centromeres_link}\n") if centromeres else None
            config.write(f"genes: {genes_link}\n") if genes else None
            config.write(f"exclude_chr: {exclude_link}\n") if exclude_chr else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"heterozygosity: {heterozygosity}\n")
        config.write(f"workflow_call: {command}\n")

    print_onstart(printmsg.rstrip("\n"), "simulate variants: translocation")
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)


simulate.add_command(linkedreads)
simulate.add_command(snpindel)
simulate.add_command(inversion)
simulate.add_command(cnv)
simulate.add_command(translocation)
