"""Harpy workflows to simulate genomic variants and linked-reads"""
import os
import sys
import yaml
from pathlib import Path
import rich_click as click
from ._cli_types_generic import HPCProfile, InputFile, SnakemakeParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake, SNAKEMAKE_CMD
from ._misc import fetch_rule, fetch_script, snakemake_log
from ._printing import print_error, workflow_info
from ._validations import check_fasta, validate_barcodefile

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
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
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
            "name": "Diploid Options",
            "options": ["--heterozygosity", "--only-vcf"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--container", "--hpc", "--output-dir", "--prefix", "--quiet", "--random-seed", "--snakemake", "--help"],
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
            "name": "Diploid Options",
            "options": ["--heterozygosity", "--only-vcf"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--container", "--hpc", "--output-dir", "--prefix", "--quiet", "--random-seed", "--snakemake", "--help"],
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
            "name": "Diploid Options",
            "options": ["--heterozygosity", "--only-vcf"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--container", "--hpc", "--output-dir", "--prefix", "--quiet", "--random-seed", "--snakemake", "--help"],
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
            "name": "Diploid Options",
            "options": ["--heterozygosity", "--only-vcf"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--container", "--hpc", "--output-dir", "--prefix", "--quiet", "--random-seed", "--snakemake", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/simulate/simulate-linkedreads")
@click.option('-b', '--barcodes', type = click.Path(exists=True, dir_okay=False, readable=True), help = "File of linked-read barcodes to add to reads")
@click.option('-s', '--distance-sd', type = click.IntRange(min = 1), default = 15, show_default=True,  help = "Standard deviation of read-pair distance")
@click.option('-m', '--molecules-per', type = click.IntRange(min = 1, max = 4700), default = 10, show_default=True,  help = "Average number of molecules per partition")
@click.option('-l', '--molecule-length', type = click.IntRange(min = 2), default = 100, show_default=True,  help = "Mean molecule length (kbp)")
@click.option('-r', '--mutation-rate', type = click.FloatRange(min = 0), default=0.001, show_default=True,  help = "Random mutation rate for simulating reads")
@click.option('-d', '--outer-distance', type = click.IntRange(min = 100), default = 350, show_default= True, help = "Outer distance between paired-end reads (bp)")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/linkedreads", help = 'Output directory name')
@click.option('-p', '--partitions', type = click.IntRange(min = 1), default=1500, show_default=True,  help = "Number (in thousands) of partitions/beads to generate")
@click.option('-n', '--read-pairs', type = click.FloatRange(min = 0.001), default = 600, show_default=True,  help = "Number (in millions) of read pairs to simulate")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome_hap1', required=True, type = InputFile("fasta", gzip_ok = True), nargs=1)
@click.argument('genome_hap2', required=True, type = InputFile("fasta", gzip_ok = True), nargs=1)
def linkedreads(genome_hap1, genome_hap2, output_dir, outer_distance, mutation_rate, distance_sd, barcodes, read_pairs, molecule_length, partitions, molecules_per, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Create linked reads from a genome
 
    Two haplotype genomes (un/compressed fasta) need to be provided as inputs at the end of the command. If
    you don't have a diploid genome, you can simulate one with `harpy simulate` as described [in the documentation](https://pdimens.github.io/harpy/blog/simulate_diploid/).

    If not providing a file for `--barcodes`, Harpy will generate a file containing the original
    (96^4) set of 24-basepair haplotagging barcodes (~2GB disk space). The `--barcodes` file is expected to have one
    linked-read barcode per line, given as nucleotides.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/simulate_linkedreads.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        command += f" --workflow-profile {hpc}"
    if snakemake:
        command += f" {snakemake}"

    check_fasta(genome_hap1)
    check_fasta(genome_hap2)
    if barcodes:
        bc_len = validate_barcodefile(barcodes, True, quiet)
    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fetch_rule(workflowdir, "simulate_linkedreads.smk")
    fetch_script(workflowdir, "HaploSim.pl")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_linkedreads")
    conda_envs = ["simulations"]
    configs = {
        "workflow" : "simulate linkedreads",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "outer_distance" : outer_distance,
        "distance_sd" : distance_sd,
        "read_pairs" : read_pairs,
        "mutation_rate" : mutation_rate,
        "molecule_length" : molecule_length,
        "partitions" : partitions,
        "molecules_per_partition" : molecules_per,
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        'barcodes': {
            "file": Path(barcodes).resolve().as_posix() if barcodes else f"{workflowdir}/input/haplotag_barcodes.txt",
            "length": bc_len if barcodes else 24
        },
        "inputs" : {
            "genome_hap1" : Path(genome_hap1).resolve().as_posix(),
            "genome_hap2" : Path(genome_hap2).resolve().as_posix(),
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Genome Haplotype 1:", os.path.basename(genome_hap1)),
        ("Genome Haplotype 2:", os.path.basename(genome_hap2)),
        ("Barcodes:", os.path.basename(barcodes) if barcodes else "Haplotagging Default"),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "simulate_linkedreads", start_text, output_dir, sm_log, quiet, "workflow/simulate.reads.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "This workflow can be quite technical, please read the docs for more information: https://pdimens.github.io/harpy/workflows/simulate/simulate-variants")
@click.option('-s', '--snp-vcf', type=InputFile("vcf", gzip_ok = False), help = 'VCF file of known snps to simulate')
@click.option('-i', '--indel-vcf', type=InputFile("vcf", gzip_ok = False), help = 'VCF file of known indels to simulate')
@click.option('-n', '--snp-count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random snps to simluate")
@click.option('-m', '--indel-count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random indels to simluate")
@click.option('-r', '--titv-ratio', type = click.FloatRange(min=0), default= 0.5, show_choices=True, show_default=True, help = "Transition/Transversion ratio for snps")
@click.option('-d', '--indel-ratio', type = click.FloatRange(min=0), default= 1, show_choices=True, show_default=True, help = "Insertion/Deletion ratio for indels")
@click.option('-a', '--indel-size-alpha', type = click.FloatRange(min=0), hidden = True, default = 2.0, help = "Exponent Alpha for power-law-fitted size distribution")
@click.option('-l', '--indel-size-constant', type = click.FloatRange(min=0), default = 0.5, hidden = True, help = "Exponent constant for power-law-fitted size distribution")
@click.option('-c', '--centromeres', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of centromeres to avoid")
@click.option('-y', '--snp-gene-constraints', type = click.Choice(["noncoding", "coding", "2d", "4d"]), help = "How to constrain randomly simulated SNPs {`noncoding`,`coding`,`2d`,`4d`}")
@click.option('-g', '--genes', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = 'heterozygosity to simulate diploid variants')
@click.option('--only-vcf',  is_flag = True, default = False, help = 'If setting heterozygosity, only create the vcf rather than the fasta files')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/snpindel", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prefix', type = str, default= "sim", show_default=True, help = "Naming prefix for output files")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome', required=True, type=InputFile("fasta", gzip_ok = True), nargs=1)
def snpindel(genome, snp_vcf, indel_vcf, only_vcf, output_dir, prefix, snp_count, indel_count, titv_ratio, indel_ratio, indel_size_alpha, indel_size_constant, centromeres, genes, snp_gene_constraints, heterozygosity, exclude_chr, random_seed, snakemake, quiet, hpc, container, setup_only):
    """
    Introduce snps and/or indels into a genome
 
    ### Haploid
    Use either a VCF file to simulate known variants via `--snp-vcf` and/or `--indel-vcf` or the command line options listed below
    to simulate random variants of the selected type.
    
    ### Diploid
    To simulate a diploid genome with heterozygous and homozygous variants, set `--heterozygosity` to a value greater than `0`.
    Use `--only-vcf` alongside `--heterozygosity` if you want to generate only the variant VCF file(s) without creating the diploid genome.
    
    The ratio parameters control different things for snp and indel variants and have special
    meanings when setting the value to either `9999` or `0` :
    
    | ratio | meaning | `9999` | `0` |
    |:---- |:---- |:--- |:---- |
    | `--snp-ratio`   | transitions / transversions | transit. only | transv. only |
    | `--indel-ratio` | insertions / deletions      | insert. only  | delet. only  |
    """
    if (snp_gene_constraints and not genes) or (genes and not snp_gene_constraints):
        print_error("missing option", "The options `--genes` and `--snp-coding-partition` must be used together for SNP variants.")
        sys.exit(1)
    if (not snp_vcf and snp_count == 0) and (not indel_vcf and indel_count == 0):
        print_error("missing option", "You must either provide a vcf file of known variants to simulate or a count of that variant to randomly simulate.")
        sys.exit(1)
    if snp_count > 0 and snp_vcf:
        print_error("conflicting arguments", "You can either simulate random SNPs (--snp-count) or known snps (--snp-vcf), but not both.")
        sys.exit(1)    
    if indel_count > 0 and indel_vcf:
        print_error("conflicting arguments", "You can either simulate random indels (--indel-count) or known indels (--indel-vcf), but not both.")
        sys.exit(1)

    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores 2'
    command += f" --snakefile {workflowdir}/simulate_snpindel.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        command += f" --workflow-profile {hpc}"
    if snakemake:
        command += f" {snakemake}"
    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    check_fasta(genome)
    fetch_rule(workflowdir, "simulate_snpindel.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_snpindel")
    conda_envs = ["simulations"]
    configs = {
        "workflow" : "simulate snpindel",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "prefix" : prefix,
        **({"random_seed" : random_seed} if random_seed else {}),
        "heterozygosity" : {
            "ratio" : heterozygosity,
            "only_vcf" : only_vcf,
        },
        "snp" : {
            **({"vcf" : Path(snp_vcf).resolve().as_posix()} if snp_vcf else {}),
            **({'count': snp_count} if snp_count and not snp_vcf else {}),
            **({"gene_constraints":  snp_gene_constraints} if snp_gene_constraints and not snp_vcf else {}),
            **({"titv_ratio" : titv_ratio} if titv_ratio and not snp_vcf else {})
        },
        "indel" : {
            **({"vcf" : Path(indel_vcf).resolve().as_posix()} if indel_vcf else {}),
            **({"count" : indel_count} if indel_count and not indel_vcf else {}),
            **({"indel_ratio" : indel_ratio} if indel_ratio and not indel_vcf else {}),
            **({"size_alpha" : indel_size_alpha} if indel_size_alpha and not indel_vcf else {}),
            **({"size_constant" : indel_size_constant} if indel_size_constant and not indel_vcf else {})
        },
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "inputs" : {
            "genome" : Path(genome).resolve().as_posix(),
            **({"centromeres" : Path(centromeres).resolve().as_posix()} if centromeres else {}),
            **({"genes" : Path(genes).resolve().as_posix()} if genes else {}),
            **({"excluded_chromosomes" : Path(exclude_chr).resolve().as_posix()} if exclude_chr else {})
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)
    start_text = workflow_info(
        ("Input Genome:", os.path.basename(genome)),
        ("SNP File:", os.path.basename(snp_vcf)) if snp_vcf else None,
        ("Random SNPs:", snp_count) if snp_count > 0 else None,
        ("Indel File:", os.path.basename(indel_vcf)) if indel_vcf else None,
        ("Random Indels:", f"{indel_count}") if indel_count else None,
        ("Centromere GFF:", os.path.basename(centromeres)) if centromeres else None,
        ("Genes GFF:", os.path.basename(genes)) if genes else None,
        ("Excluded Chromosomes:", os.path.basename(exclude_chr)) if exclude_chr else None,
        ("Heterozygosity:", heterozygosity) if heterozygosity > 0 else None,
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "simulate_snpindel", start_text, output_dir, sm_log, quiet, "workflow/simulate.snpindel.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Please Documentation: https://pdimens.github.io/harpy/workflows/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True), help = 'VCF file of known inversions to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random inversions to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum inversion size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum inversion size (bp)")
@click.option('-c', '--centromeres', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = 'heterozygosity to simulate diploid variants')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True), help = "Text file of chromosomes to avoid")
@click.option('-p', '--prefix', type = str, default= "sim", show_default=True, help = "Naming prefix for output files")
@click.option('--only-vcf',  is_flag = True, default = False, help = 'If setting heterozygosity, only create the vcf rather than the fasta files')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/inversion", show_default=True,  help = 'Output directory name')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome', required=True, type=InputFile("fasta", gzip_ok = True), nargs=1)
def inversion(genome, vcf, only_vcf, prefix, output_dir, count, min_size, max_size, centromeres, genes, heterozygosity, exclude_chr, random_seed, snakemake, quiet, hpc, container, setup_only):
    """
    Introduce inversions into a genome
 
    ### Haploid
    Use either a VCF file to simulate known inversions or the command line options listed below to simulate random inversions.

    ### Diploid
    To simulate a diploid genome with heterozygous and homozygous variants, set `--heterozygosity` to a value greater than `0`.
    Use `--only-vcf` alongside `--heterozygosity` if you want to generate only the variant VCF file(s) without creating the diploid genome.
    """
    if not vcf and count == 0:
        print_error("missing option", "Provide either a `--count` of cnv to randomly simulate or a `--vcf` of known variants to simulate.")
        sys.exit(1)
    if vcf and count > 0:
        print_error("conflicting arguments", "You can either simulate random inversions (--count) or known inversions (--vcf), but not both.")
        sys.exit(1)   

    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores 2'
    command += f" --snakefile {workflowdir}/simulate_variants.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        command += f" --workflow-profile {hpc}"
    if snakemake:
        command += f" {snakemake}"

    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    check_fasta(genome)
    fetch_rule(workflowdir, "simulate_variants.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_inversion")
    conda_envs = ["simulations"]
    configs = {
        "workflow" : "simulate inversion",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "prefix" : prefix,
        **({"random_seed" : random_seed} if random_seed else {}),
        "heterozygosity" : {
            "ratio" : heterozygosity,
            "only_vcf" : only_vcf,
        },
        "inversion" : {
            **({"vcf" : Path(vcf).resolve().as_posix()} if vcf else {}),
            **({'count': count} if not vcf else {}),
            **({"min_size":  min_size} if not vcf else {}),
            **({"max_size" : max_size} if not vcf else {})
        },
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "inputs" : {
            "genome" : Path(genome).resolve().as_posix(),
            **({"centromeres" : Path(centromeres).resolve().as_posix()} if centromeres else {}),
            **({"genes" : Path(genes).resolve().as_posix()} if genes else {}),
            **({"excluded_chromosomes" : Path(exclude_chr).resolve().as_posix()} if exclude_chr else {})
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)
    start_text = workflow_info(
        ("Input Genome:", os.path.basename(genome)),
        ("Inversion File:", os.path.basename(vcf)) if vcf else ("Random Inversions:", count),
        ("Centromere GFF:", os.path.basename(centromeres)) if centromeres else None,
        ("Genes GFF:", os.path.basename(genes)) if genes else None,
        ("Excluded Chromosomes:", os.path.basename(exclude_chr)) if exclude_chr else None,
        ("Heterozygosity:", heterozygosity) if heterozygosity > 0 else None,
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "simulate_inversion", start_text, output_dir, sm_log, quiet, "workflow/simulate.inversion.summary")


@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Please Documentation: https://pdimens.github.io/harpy/workflows/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True), help = 'VCF file of known copy number variants to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random variants to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum variant size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum variant size (bp)")
@click.option('-d', '--dup-ratio', type = click.FloatRange(min=0), default = 1, show_choices=True, show_default=True, help = "Ratio for the selected variant")
@click.option('-l', '--gain-ratio', type = click.FloatRange(min=0), default=1, show_default=True, help = " Relative ratio of DNA gain over DNA loss")
@click.option('-y', '--max-copy', type = click.IntRange(min = 1), default=10, show_default=True, help = "Maximum number of copies")
@click.option('-c', '--centromeres', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = 'heterozygosity to simulate diploid variants')
@click.option('--only-vcf',  is_flag = True, default = False, help = 'If setting heterozygosity, only create the vcf rather than the fasta files')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/cnv", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prefix', type = str, default= "sim", show_default=True, help = "Naming prefix for output files")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome', required=True, type=InputFile("fasta", gzip_ok = True), nargs=1)
def cnv(genome, output_dir, vcf, only_vcf, prefix, count, min_size, max_size, dup_ratio, max_copy, gain_ratio, centromeres, genes, heterozygosity, exclude_chr, random_seed, snakemake, quiet, hpc, container, setup_only):
    """
    Introduce copy number variants into a genome
 
    ### Haploid
    Use either a VCF file to simulate known CNVs or the command line options listed below to simulate random CNVs.
      
    ### Diploid
    To simulate a diploid genome with heterozygous and homozygous variants, set `--heterozygosity` to a value greater than `0`.
    Use `--only-vcf` alongside `--heterozygosity` if you want to generate only the variant VCF file(s) without creating the diploid genome.
 
    The two ratio parameters control different things and have special meanings when setting their value to either `9999` or `0`:
    
    | ratio | meaning | `9999` | `0` |
    |:----|:----|:---|:----|
    | `--dup-ratio` | tandem / dispersed | tand. only | disp. only |
    | `--gain-ratio` | copy gain / loss | gain only | loss only| 
    """
    if not vcf and count == 0:
        print_error("missing option", "Provide either a `--count` of cnv to randomly simulate or a `--vcf` of known cnv to simulate.")
        sys.exit(1)
    if vcf and count > 0:
        print_error("conflicting arguments", "You can either simulate random CNVs (--count) or known CNVs (--vcf), but not both.")
        sys.exit(1) 
    
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores 2'
    command += f" --snakefile {workflowdir}/simulate_variants.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        command += f" --workflow-profile {hpc}"
    if snakemake:
        command += f" {snakemake}"

    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    check_fasta(genome)
    fetch_rule(workflowdir, "simulate_variants.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_cnv")
    conda_envs = ["simulations"]
    configs = {
        "workflow" : "simulate cnv",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "prefix" : prefix,
        **({"random_seed" : random_seed} if random_seed else {}),
        "heterozygosity" : {
            "ratio" : heterozygosity,
            "only_vcf" : only_vcf,
        },
        "cnv" : {
            **({"vcf" : Path(vcf).resolve().as_posix()} if vcf else {}),
            **({'count': count} if not vcf else {}),
            **({"min_size":  min_size} if not vcf else {}),
            **({"max_size" : max_size} if not vcf else {}),
            **({"duplication_ratio" : dup_ratio} if not vcf else {}),
            **({"max_copy" : max_copy} if not vcf else {}),
            **({"gain_ratio" : gain_ratio} if not vcf else {})
        },
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "inputs" : {
            "genome" : Path(genome).resolve().as_posix(),
            **({"centromeres" : Path(centromeres).resolve().as_posix()} if centromeres else {}),
            **({"genes" : Path(genes).resolve().as_posix()} if genes else {}),
            **({"excluded_chromosomes" : Path(exclude_chr).resolve().as_posix()} if exclude_chr else {})
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)
    start_text = workflow_info(
        ("Input Genome:", os.path.basename(genome)),
        ("CNV File:", os.path.basename(vcf)) if vcf else ("Random CNVs:", count),
        ("Centromere GFF:", os.path.basename(centromeres)) if centromeres else None,
        ("Genes GFF:", os.path.basename(genes)) if genes else None,
        ("Excluded Chromosomes:", os.path.basename(exclude_chr)) if exclude_chr else None,
        ("Heterozygosity:", heterozygosity) if heterozygosity > 0 else None,
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "simulate_cnv", start_text, output_dir, sm_log, quiet, "workflow/simulate.cnv.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Please Documentation: https://pdimens.github.io/harpy/workflows/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True), help = 'VCF file of known translocations to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random translocations to simluate")
@click.option('-c', '--centromeres', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = 'heterozygosity to simulate diploid variants')
@click.option('--only-vcf',  is_flag = True, default = False, help = 'If setting heterozygosity, only create the vcf rather than the fasta files')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/translocation", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prefix', type = str, default= "sim", show_default=True, help = "Naming prefix for output files")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for simulation")
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome', required=True, type=InputFile("fasta", gzip_ok = True), nargs=1)
def translocation(genome, output_dir, prefix, vcf, only_vcf, count, centromeres, genes, heterozygosity, exclude_chr, random_seed, snakemake, quiet, hpc, container, setup_only):
    """
    Introduce translocations into a genome
 
    ### Haploid
    Use either a VCF file to simulate known translocations or the command line options listed below to simulate random translocations.
      
    ### Diploid
    To simulate a diploid genome with heterozygous and homozygous variants, set `--heterozygosity` to a value greater than `0`.
    Use `--only-vcf` alongside `--heterozygosity` if you want to generate only the variant VCF file(s) without creating the diploid genome.
    """
    if not vcf and count == 0:
        print_error("missing option", "Provide either a `--count` of cnv to randomly simulate or a `--vcf` of known cnv to simulate.")
        sys.exit(1)
    if vcf and count > 0:
        print_error("conflicting arguments", "You can either simulate random translocations (--count) or known translocations (--vcf), but not both.")
        sys.exit(1) 
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores 2'
    command += f" --snakefile {workflowdir}/simulate_variants.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        command += f" --workflow-profile {hpc}"
    if snakemake:
        command += f" {snakemake}"

    # instantiate workflow directory
    # move necessary files to workflow dir
    os.makedirs(f"{workflowdir}/input/", exist_ok= True)
    check_fasta(genome)
    fetch_rule(workflowdir, "simulate_variants.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_translocation")
    conda_envs = ["simulations"]
    configs = {
        "workflow" : "simulate translocation",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "prefix" : prefix,
        **({"random_seed" : random_seed} if random_seed else {}),
        "heterozygosity" : {
            "ratio" : heterozygosity,
            "only_vcf" : only_vcf
        },
        "translocation" : {
            **({"vcf" : Path(vcf).resolve().as_posix()} if vcf else {}),
            **({'count': count} if not vcf else {})
        },
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "inputs" : {
            "genome" : Path(genome).resolve().as_posix(),
            **({"centromeres" : Path(centromeres).resolve().as_posix()} if centromeres else {}),
            **({"genes" : Path(genes).resolve().as_posix()} if genes else {}),
            **({"excluded_chromosomes" : Path(exclude_chr).resolve().as_posix()} if exclude_chr else {})
        }
    }
    # setup the config file depending on inputs
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)
    start_text = workflow_info(
        ("Input Genome:", os.path.basename(genome)),
        ("Translocation File:", os.path.basename(vcf)) if vcf else ("Random Translocations:", f"{count}"),
        ("Centromere GFF:", os.path.basename(centromeres)) if centromeres else None,
        ("Genes GFF:", os.path.basename(genes)) if genes else None,
        ("Excluded Chromosomes:", os.path.basename(exclude_chr)) if exclude_chr else None,
        ("Heterozygosity:", heterozygosity) if heterozygosity > 0 else None,
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "simulate_translocation", start_text, output_dir, sm_log, quiet, "workflow/simulate.translocation.summary")

simulate.add_command(linkedreads)
simulate.add_command(snpindel)
simulate.add_command(inversion)
simulate.add_command(cnv)
simulate.add_command(translocation)
