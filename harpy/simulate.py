"""Harpy workflows to simulate genomic variants and linked-reads"""
import os
import sys
from pathlib import Path
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, fetch_script, snakemake_log, launch_snakemake
from .printfunctions import print_error
from .validations import validate_input_by_ext, check_fasta

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def simulate():
    """
    Simulate variants or linked-reads from a genome

    To simulate genomic variants, provide an additional subcommand {`snpindel`,`inversion`,`cnv`,`translocation`} 
    to get more information about that workflow. The variant simulator (`simuG`) can only simulate
    one type of variant at a time, so you may need to run it a few times if you want multiple variant types.
    Use `simulate linkedreads` to simulate haplotag linked-reads from a diploid genome, which you can create by simulating
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
            "commands": ["cnv", "inversion", "snpindel", "translocation"],
        }
    ]
}

docstring = {
    "harpy simulate linkedreads": [
        {
            "name": "Parameters",
            "options": ["--barcodes", "--distance-sd", "--outer-distance", "--molecule-length", "--molecules-per", "--mutation-rate", "--partitions", "--read-pairs"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
        },     
    ],
    "harpy simulate snpindel": [
        {
            "name": "Known Variants",
            "options": ["--indel-vcf", "--snp-vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--centromeres", "--exclude-chr", "--genes", "--indel-count", "--indel-ratio", "--snp-count", "--snp-gene-constraints", "--titv-ratio"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--heterozygosity", "--hpc", "--output-dir", "--prefix", "--quiet", "--randomseed", "--snakemake", "--help"],
        },
    ],
    "harpy simulate inversion": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--centromeres", "--count", "--exclude-chr", "--genes", "--max-size", "--min-size"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--heterozygosity", "--hpc", "--output-dir", "--prefix", "--quiet", "--randomseed", "--snakemake", "--help"],
        },
    ],
    "harpy simulate cnv": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--centromeres", "--count", "--dup-ratio", "--exclude-chr", "--gain-ratio", "--genes",  "--max-copy", "--max-size", "--min-size"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--heterozygosity", "--hpc", "--output-dir", "--prefix", "--quiet", "--randomseed", "--snakemake", "--help"],
        },
    ],
    "harpy simulate translocation": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--centromeres", "--count", "--exclude-chr", "--genes"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--heterozygosity", "--hpc", "--output-dir", "--prefix", "--quiet", "--randomseed", "--snakemake", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-linkedreads")
@click.option('-b', '--barcodes', type = click.Path(exists=True, dir_okay=False), help = "File of linked-read barcodes to add to reads")
@click.option('-s', '--distance-sd', type = click.IntRange(min = 1), default = 15, show_default=True,  help = "Standard deviation of read-pair distance")
@click.option('-m', '--molecules-per', type = click.IntRange(min = 1), default = 10, show_default=True,  help = "Average number of molecules per partition")
@click.option('-l', '--molecule-length', type = click.IntRange(min = 2), default = 100, show_default=True,  help = "Mean molecule length (kbp)")
@click.option('-r', '--mutation-rate', type = click.FloatRange(min = 0), default=0.001, show_default=True,  help = "Random mutation rate for simulating reads")
@click.option('-d', '--outer-distance', type = click.IntRange(min = 100), default = 350, show_default= True, help = "Outer distance between paired-end reads (bp)")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/linkedreads", help = 'Output directory name')
@click.option('-p', '--partitions', type = click.IntRange(min = 1), default=1500, show_default=True,  help = "Number (in thousands) of partitions/beads to generate")
@click.option('-n', '--read-pairs', type = click.FloatRange(min = 0.001), default = 600, show_default=True,  help = "Number (in millions) of read pairs to simulate")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--config-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Create the config.yaml file and exit')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome_hap1', required=True, type=click.Path(exists=True, dir_okay=False, readable=True), nargs=1)
@click.argument('genome_hap2', required=True, type=click.Path(exists=True, dir_okay=False, readable=True), nargs=1)
def linkedreads(genome_hap1, genome_hap2, output_dir, outer_distance, mutation_rate, distance_sd, barcodes, read_pairs, molecule_length, partitions, molecules_per, threads, snakemake, quiet, hpc, conda, config_only):
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
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock  --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/simulate_linkedreads.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    check_fasta(genome_hap1)
    check_fasta(genome_hap2)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fetch_rule(workflowdir, "simulate_linkedreads.smk")
    fetch_script(workflowdir, "LRSIM_harpy.pl")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_linkedreads")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: simulate linkedreads\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        if barcodes:
            config.write(f"barcodes: {Path(barcodes).resolve()}\n")
        config.write(f"outer_distance: {outer_distance}\n")
        config.write(f"distance_sd: {distance_sd}\n")
        config.write(f"read_pairs: {read_pairs}\n")
        config.write(f"mutation_rate: {mutation_rate}\n")
        config.write(f"molecule_length: {molecule_length}\n")
        config.write(f"partitions: {partitions}\n")
        config.write(f"molecules_per_partition: {molecules_per}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome_hap1: {Path(genome_hap1).resolve()}\n")
        config.write(f"  genome_hap2: {Path(genome_hap2).resolve()}\n")
    if config_only:
        sys.exit(0)

    generate_conda_deps()
    start_text = f"Genome Haplotype 1: {os.path.basename(genome_hap1)}\n"
    start_text += f"Genome Haplotype 2: {os.path.basename(genome_hap2)}\n"
    start_text += f"Barcodes: {os.path.basename(barcodes)}\n" if barcodes else "Barcodes: 10X Default\n"
    start_text += f"Output Directory: {output_dir}/"
    launch_snakemake(command, "simulate_linkedreads", start_text, output_dir, sm_log)

@click.command(no_args_is_help = True, epilog = "This workflow can be quite technical, please read the docs for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-variants")
@click.option('-s', '--snp-vcf', type=click.Path(exists=True, dir_okay=False, readable=True), help = 'VCF file of known snps to simulate')
@click.option('-i', '--indel-vcf', type=click.Path(exists=True, dir_okay=False, readable=True), help = 'VCF file of known indels to simulate')
@click.option('-n', '--snp-count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random snps to simluate")
@click.option('-m', '--indel-count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random indels to simluate")
@click.option('-r', '--titv-ratio', type = click.FloatRange(min=0), default= 0.5, show_choices=True, show_default=True, help = "Transition/Transversion ratio for snps")
@click.option('-d', '--indel-ratio', type = click.FloatRange(min=0), default= 1, show_choices=True, show_default=True, help = "Insertion/Deletion ratio for indels")
@click.option('-a', '--indel-size-alpha', type = click.FloatRange(min=0), hidden = True, default = 2.0, help = "Exponent Alpha for power-law-fitted size distribution")
@click.option('-l', '--indel-size-constant', type = click.FloatRange(min=0), default = 0.5, hidden = True, help = "Exponent constant for power-law-fitted size distribution")
@click.option('-c', '--centromeres', type = click.Path(exists=True, dir_okay=False, readable=True), help = "GFF3 file of centromeres to avoid")
@click.option('-y', '--snp-gene-constraints', type = click.Choice(["noncoding", "coding", "2d", "4d"]), help = "How to constrain randomly simulated SNPs {`noncoding`,`coding`,`2d`,`4d`}")
@click.option('-g', '--genes', type = click.Path(exists=True, readable=True), help = "GFF3 file of genes to use with `--snp-gene-constraints`")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/snpindel", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prefix', type = str, default= "sim.snpindel", show_default=True, help = "Naming prefix for output files")
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--config-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Create the config.yaml file and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome', required=True, type=click.Path(exists=True, dir_okay=False, readable=True), nargs=1)
def snpindel(genome, snp_vcf, indel_vcf, output_dir, prefix, snp_count, indel_count, titv_ratio, indel_ratio, indel_size_alpha, indel_size_constant, centromeres, genes, snp_gene_constraints, heterozygosity, exclude_chr, randomseed, snakemake, quiet, hpc, conda, config_only):
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
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores 1 --directory . '
    command += f"--snakefile {workflowdir}/simulate_snpindel.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)   
    check_fasta(genome)
    start_text = f"Inpute Genome: {os.path.basename(genome)}\nOutput Directory: {output_dir}/\n"
    if snp_vcf:
        validate_input_by_ext(snp_vcf, "--snp-vcf", ["vcf","vcf.gz","bcf"])
        start_text += f"SNPs: from vcf ({os.path.basename(snp_vcf)})\n"
    elif snp_count > 0:
        start_text += "SNPs: random\n"
    if indel_vcf:
        validate_input_by_ext(indel_vcf, "--indel-vcf", ["vcf","vcf.gz","bcf"])
        start_text += f"Indels: from vcf: ({os.path.basename(indel_vcf)})\n"
    elif indel_count > 0:
        start_text += "Indels: random\n"
    if centromeres:
        validate_input_by_ext(centromeres, "--centromeres", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        start_text += f"Centromere GFF: {os.path.basename(centromeres)}\n"
    if genes:
        validate_input_by_ext(genes, "--genes", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        start_text += f"Genes GFF: {os.path.basename(genes)}\n"
    if exclude_chr:
        start_text += f"Excluded Chromosomes: {os.path.basename(exclude_chr)}\n"
    fetch_rule(workflowdir, "simulate_snpindel.smk")
    fetch_script(workflowdir, "simuG.pl")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_snpindel")

    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: simulate snpindel\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"prefix: {prefix}\n")
        config.write(f"heterozygosity: {heterozygosity}\n")
        if not snp_vcf:
            config.write(f"snp_count: {snp_count}\n") if snp_count else None
            config.write(f"snp_gene_constraints: {snp_gene_constraints}\n") if snp_gene_constraints else None
            config.write(f"titv_ratio: {titv_ratio}\n") if titv_ratio else None
        if not indel_vcf:
            config.write(f"indel_count: {snp_count}\n") if indel_count else None
            config.write(f"indel_ratio: {indel_ratio}\n") if indel_ratio else None
            config.write(f"indel_size_alpha: {indel_size_alpha}\n") if indel_size_alpha else None
            config.write(f"indel_size_constant: {indel_size_constant}\n") if indel_size_constant else None
        config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome: {Path(genome).resolve()}\n")
        if snp_vcf:
            config.write(f"  snp_vcf: {Path(snp_vcf).resolve()}\n")
        if indel_vcf:
            config.write(f"  indel_vcf: {Path(indel_vcf).resolve()}\n")
        config.write(f"  centromeres: {Path(centromeres).resolve()}\n") if centromeres else None
        config.write(f"  genes: {Path(genes).resolve()}\n") if genes else None
        config.write(f"  exclude_chr: {Path(exclude_chr).resolve()}\n") if exclude_chr else None
    if config_only:
        sys.exit(0)

    generate_conda_deps()
    launch_snakemake(command, "simulate_snpindel", start_text.rstrip("\n"), output_dir, sm_log)


@click.command(no_args_is_help = True, epilog = "Please See the documentation for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True), help = 'VCF file of known inversions to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random inversions to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum inversion size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum inversion size (bp)")
@click.option('-c', '--centromeres', type = click.Path(exists=True, dir_okay=False, readable=True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True, dir_okay=False, readable=True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True), help = "Text file of chromosomes to avoid")
@click.option('-p', '--prefix', type = str, default= "sim.inversion", show_default=True, help = "Naming prefix for output files")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/inversion", show_default=True,  help = 'Output directory name')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--config-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Create the config.yaml file and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome', required=True, type=click.Path(exists=True, dir_okay=False, readable=True), nargs=1)
def inversion(genome, vcf, prefix, output_dir, count, min_size, max_size, centromeres, genes, heterozygosity, exclude_chr, randomseed, snakemake, quiet, hpc, conda, config_only):
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
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores 1 --directory . '
    command += f"--snakefile {workflowdir}/simulate_variants.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    check_fasta(genome)
    start_text = f"Inpute Genome: {os.path.basename(genome)}\nOutput Directory: {output_dir}/\n"

    if vcf:
        validate_input_by_ext(vcf, "--vcf", ["vcf","vcf.gz","bcf"])
        start_text += f"Input VCF: {os.path.basename(vcf)}\n"
    else:
        start_text += "Mode: Random variants\n"
    if centromeres:
        validate_input_by_ext(centromeres, "--centromeres", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        start_text += f"Centromere GFF: {os.path.basename(centromeres)}\n"
    if genes:
        validate_input_by_ext(genes, "--genes", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        start_text += f"Genes GFF: {os.path.basename(genes)}\n"
    if exclude_chr:
        start_text += f"Excluded Chromosomes: {os.path.basename(exclude_chr)}\n"
    fetch_rule(workflowdir, "simulate_variants.smk")
    fetch_script(workflowdir, "simuG.pl")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_inversion")

    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: simulate inversion\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write("variant_type: inversion\n")
        config.write(f"prefix: {prefix}\n")
        config.write(f"heterozygosity: {heterozygosity}\n")
        if not vcf:
            config.write(f"count: {count}\n")
            config.write(f"min_size: {min_size}\n") if min_size else None
            config.write(f"max_size: {max_size}\n") if max_size else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome: {Path(genome).resolve()}\n")
        if vcf:
            config.write(f"  vcf: {Path(vcf).resolve()}\n")
        else:
            config.write(f"  centromeres: {Path(centromeres).resolve()}\n") if centromeres else None
            config.write(f"  genes: {Path(genes).resolve()}\n") if genes else None
            config.write(f"  exclude_chr: {Path(exclude_chr).resolve()}\n") if exclude_chr else None
    if config_only:
        sys.exit(0)

    generate_conda_deps()
    launch_snakemake(command, "simulate_inversion", start_text.rstrip("\n"), output_dir, sm_log)


@click.command(no_args_is_help = True, epilog = "Please See the documentation for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True), help = 'VCF file of known copy number variants to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random variants to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum variant size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum variant size (bp)")
@click.option('-d', '--dup-ratio', type = click.FloatRange(min=0), default = 1, show_choices=True, show_default=True, help = "Ratio for the selected variant")
@click.option('-l', '--gain-ratio', type = click.FloatRange(min=0), default=1, show_default=True, help = " Relative ratio of DNA gain over DNA loss")
@click.option('-y', '--max-copy', type = click.IntRange(min = 1), default=10, show_default=True, help = "Maximum number of copies")
@click.option('-c', '--centromeres', type = click.Path(exists=True, dir_okay=False, readable=True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True, dir_okay=False, readable=True), help = "GFF3 file of genes to avoid when simulating (requires `--snp-coding-partition` for SNPs)")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/cnv", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prefix', type = str, default= "sim.cnv", show_default=True, help = "Naming prefix for output files")
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--config-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Create the config.yaml file and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome', required=True, type=click.Path(exists=True, dir_okay=False, readable=True), nargs=1)
def cnv(genome, output_dir, vcf, prefix, count, min_size, max_size, dup_ratio, max_copy, gain_ratio, centromeres, genes, heterozygosity, exclude_chr, randomseed, snakemake, quiet, hpc, conda, config_only):
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
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores 1 --directory . '
    command += f"--snakefile {workflowdir}/simulate_variants.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    check_fasta(genome)
    start_text = f"Inpute Genome: {os.path.basename(genome)}\nOutput Directory: {output_dir}/\n"

    if vcf:
        validate_input_by_ext(vcf, "--vcf", ["vcf","vcf.gz","bcf"])
        start_text += f"Input VCF: {os.path.basename(vcf)}\n"
    else:
        start_text += "Mode: Random variants\n"
    if centromeres:
        validate_input_by_ext(centromeres, "--centromeres", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        start_text += f"Centromere GFF: {os.path.basename(centromeres)}\n"
    if genes:
        validate_input_by_ext(genes, "--genes", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        start_text += f"Genes GFF: {os.path.basename(genes)}\n"
    if exclude_chr:
        start_text += f"Excluded Chromosomes: {os.path.basename(exclude_chr)}\n"
    fetch_rule(workflowdir, "simulate_variants.smk")
    fetch_script(workflowdir, "simuG.pl")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_cnv")

    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: simulate cnv\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write("variant_type: cnv\n")
        config.write(f"prefix: {prefix}\n")
        config.write(f"heterozygosity: {heterozygosity}\n")
        if not vcf:
            config.write(f"count: {count}\n")
            config.write(f"min_size: {min_size}\n") if min_size else None
            config.write(f"max_size: {max_size}\n") if max_size else None
            config.write(f"dup_ratio: {dup_ratio}\n") if dup_ratio else None
            config.write(f"cnv_max_copy: {max_copy}\n") if max_copy else None
            config.write(f"gain_ratio: {gain_ratio}\n") if gain_ratio else None
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome: {Path(genome).resolve()}\n")
        if vcf:
            config.write(f"  vcf: {Path(vcf).resolve()}\n")
        else:
            config.write(f"  centromeres: {Path(centromeres).resolve()}\n") if centromeres else None
            config.write(f"  genes: {Path(genes).resolve()}\n") if genes else None
            config.write(f"  exclude_chr: {Path(exclude_chr).resolve()}\n") if exclude_chr else None
    if config_only:
        sys.exit(0)

    generate_conda_deps()
    launch_snakemake(command, "simulate_cnv", start_text.rstrip("\n"), output_dir, sm_log)

@click.command(no_args_is_help = True, epilog = "Please See the documentation for more information: https://pdimens.github.io/harpy/modules/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True), help = 'VCF file of known translocations to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random translocations to simluate")
@click.option('-c', '--centromeres', type = click.Path(exists=True, dir_okay=False, readable=True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = click.Path(exists=True, dir_okay=False, readable=True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = '% heterozygosity to simulate diploid later')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/translocation", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prefix', type = str, default= "sim.translocation", show_default=True, help = "Naming prefix for output files")
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--config-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Create the config.yaml file and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--randomseed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome', required=True, type=click.Path(exists=True, dir_okay=False, readable=True), nargs=1)
def translocation(genome, output_dir, prefix, vcf, count, centromeres, genes, heterozygosity, exclude_chr, randomseed, snakemake, quiet, hpc, conda, config_only):
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
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores 1 --directory . '
    command += f"--snakefile {workflowdir}/simulate_variants.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    check_fasta(genome)
    start_text = f"Inpute Genome: {os.path.basename(genome)}\nOutput Directory: {output_dir}/\n"
    
    if vcf:
        validate_input_by_ext(vcf, "--vcf", ["vcf","vcf.gz","bcf"])
        start_text += f"Input VCF: {os.path.basename(vcf)}\n"
    else:
        start_text += "Mode: Random variants\n"
    if centromeres:
        validate_input_by_ext(centromeres, "--centromeres", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        start_text += f"Centromere GFF: {os.path.basename(centromeres)}\n"
    if genes:
        validate_input_by_ext(genes, "--genes", [".gff",".gff3",".gff.gz", ".gff3.gz"])
        start_text += f"Genes GFF: {os.path.basename(genes)}\n"
    if exclude_chr:
        start_text += f"Excluded Chromosomes: {os.path.basename(exclude_chr)}\n"

    fetch_rule(workflowdir, "simulate_variants.smk")
    fetch_script(workflowdir, "simuG.pl")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_translocation")

    # setup the config file depending on inputs
    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: simulate translocation\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write("variant_type: translocation\n")
        config.write(f"prefix: {prefix}\n")
        if not vcf:
            config.write(f"count: {count}\n")
            config.write(f"randomseed: {randomseed}\n") if randomseed else None
        config.write(f"heterozygosity: {heterozygosity}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome: {Path(genome).resolve()}\n")
        if vcf:
            config.write(f"  vcf: {Path(vcf).resolve()}\n")
        else:
            config.write(f"  centromeres: {Path(centromeres).resolve()}\n") if centromeres else None
            config.write(f"  genes: {Path(genes).resolve()}\n") if genes else None
            config.write(f"  exclude_chr: {Path(exclude_chr).resolve()}\n") if exclude_chr else None
    if config_only:
        sys.exit(0)

    generate_conda_deps()
    launch_snakemake(command, "simulate_translocation", start_text.rstrip("\n"), output_dir, sm_log)

simulate.add_command(linkedreads)
simulate.add_command(snpindel)
simulate.add_command(inversion)
simulate.add_command(cnv)
simulate.add_command(translocation)
