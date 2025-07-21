"""Harpy workflows to simulate genomic variants and linked reads"""
import os
import sys
import yaml
import shutil
import rich_click as click
from ._cli_types_generic import HPCProfile, InputFile, SnakemakeParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, instantiate_dir, setup_snakemake, write_workflow_config
from ._printing import print_error, workflow_info
from ._validations import check_fasta

commandstring = {
    "harpy simulate": [
        {
            "name": "Genomic Variants",
            "commands": ["cnv", "inversion", "snpindel", "translocation"]
        }
    ]
}

docstring = {
    "harpy simulate snpindel": [
        {
            "name": "Known Variants",
            "options": ["--indel-vcf", "--snp-vcf"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Random Variants",
            "options": ["--centromeres", "--exclude-chr", "--genes", "--indel-count", "--indel-ratio", "--snp-count", "--snp-gene-constraints", "--titv-ratio"],
            "panel_styles": {"border_style": "green"}
        },
        {
            "name": "Diploid Options",
            "options": ["--heterozygosity", "--only-vcf"],
            "panel_styles": {"border_style": "yellow"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--prefix", "--quiet", "--random-seed", "--snakemake", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy simulate inversion": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Random Variants",
            "options": ["--centromeres", "--count", "--exclude-chr", "--genes", "--max-size", "--min-size"],
            "panel_styles": {"border_style": "green"}
        },
        {
            "name": "Diploid Options",
            "options": ["--heterozygosity", "--only-vcf"],
            "panel_styles": {"border_style": "yellow"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--prefix", "--quiet", "--random-seed", "--snakemake", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy simulate cnv": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Random Variants",
            "options": ["--centromeres", "--count", "--dup-ratio", "--exclude-chr", "--gain-ratio", "--genes",  "--max-copy", "--max-size", "--min-size"],
            "panel_styles": {"border_style": "green"}
        },
        {
            "name": "Diploid Options",
            "options": ["--heterozygosity", "--only-vcf"],
            "panel_styles": {"border_style": "yellow"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--prefix", "--quiet", "--random-seed", "--snakemake", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy simulate translocation": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Random Variants",
            "options": ["--centromeres", "--count", "--exclude-chr", "--genes"],
            "panel_styles": {"border_style": "green"}
        },
        {
            "name": "Diploid Options",
            "options": ["--heterozygosity", "--only-vcf"],
            "panel_styles": {"border_style": "yellow"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--prefix", "--quiet", "--random-seed", "--snakemake", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

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
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Simulate/snpindel", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prefix', type = str, default= "sim", show_default=True, help = "Naming prefix for output files")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
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
    Use `--only-vcf` alongside `--heterozygosity` to only generate the second VCF file and not simulate a second FASTA file.
    
    The ratio parameters control different things for snp and indel variants and have special
    meanings when setting the value to either `9999` or `0` :
    
    | ratio | meaning | `9999` | `0` |
    |:---- |:---- |:--- |:---- |
    | `--snp-ratio`   | transitions / transversions | transit. only | transv. only |
    | `--indel-ratio` | insertions / deletions      | insert. only  | delet. only  |
    """
    workflow = "simulate_snpindel"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow, True)
    ## checks and validations ##
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
    check_fasta(genome)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        "conda" if not container else "conda apptainer",
        output_dir,
        2,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "simulate_snpindel.smk")

    conda_envs = ["simulations"]
    configs = {
        "workflow" : workflow,
        "prefix" : prefix,
        **({"random_seed" : random_seed} if random_seed else {}),
        "heterozygosity" : {
            "ratio" : heterozygosity,
            "only_vcf" : only_vcf,
        },
        "snp" : {
            **({"vcf" : snp_vcf} if snp_vcf else {}),
            **({'count': snp_count} if snp_count and not snp_vcf else {}),
            **({"gene_constraints":  snp_gene_constraints} if snp_gene_constraints and not snp_vcf else {}),
            **({"titv_ratio" : titv_ratio} if titv_ratio and not snp_vcf else {})
        },
        "indel" : {
            **({"vcf" : indel_vcf} if indel_vcf else {}),
            **({"count" : indel_count} if indel_count and not indel_vcf else {}),
            **({"indel_ratio" : indel_ratio} if indel_ratio and not indel_vcf else {}),
            **({"size_alpha" : indel_size_alpha} if indel_size_alpha and not indel_vcf else {}),
            **({"size_constant" : indel_size_constant} if indel_size_constant and not indel_vcf else {})
        },
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "inputs" : {
            "genome" : genome,
            **({"centromeres" : centromeres} if centromeres else {}),
            **({"genes" : genes} if genes else {}),
            **({"excluded_chromosomes" : exclude_chr} if exclude_chr else {})
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Input Genome:", os.path.basename(genome)),
        ("SNP File:", os.path.basename(snp_vcf)) if snp_vcf else None,
        ("Random SNPs:", snp_count) if snp_count > 0 else None,
        ("Indel File:", os.path.basename(indel_vcf)) if indel_vcf else None,
        ("Random Indels:", indel_count) if indel_count else None,
        ("Centromere GFF:", os.path.basename(centromeres)) if centromeres else None,
        ("Genes GFF:", os.path.basename(genes)) if genes else None,
        ("Excluded Chromosomes:", os.path.basename(exclude_chr)) if exclude_chr else None,
        ("Heterozygosity:", heterozygosity) if heterozygosity > 0 else None,
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/simulate.snpindel.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Please Documentation: https://pdimens.github.io/harpy/workflows/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), help = 'VCF file of known inversions to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random inversions to simluate")
@click.option('-m', '--min-size', type = click.IntRange(min = 1), default = 1000, show_default= True, help = "Minimum inversion size (bp)")
@click.option('-x', '--max-size', type = click.IntRange(min = 1), default = 100000, show_default= True, help = "Maximum inversion size (bp)")
@click.option('-c', '--centromeres', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = 'heterozygosity to simulate diploid variants')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), help = "Text file of chromosomes to avoid")
@click.option('-p', '--prefix', type = str, default= "sim", show_default=True, help = "Naming prefix for output files")
@click.option('--only-vcf',  is_flag = True, default = False, help = 'If setting heterozygosity, only create the vcf rather than the fasta files')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Simulate/inversion", show_default=True,  help = 'Output directory name')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
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
    Use `--only-vcf` alongside `--heterozygosity` to only generate the second VCF file and not simulate a second FASTA file.
    """
    workflow = "simulate_inversion"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow, True)
    ## checks and validations ##
    if not vcf and count == 0:
        print_error("missing option", "Provide either a `--count` of cnv to randomly simulate or a `--vcf` of known variants to simulate.")
        sys.exit(1)
    if vcf and count > 0:
        print_error("conflicting arguments", "You can either simulate random inversions (--count) or known inversions (--vcf), but not both.")
        sys.exit(1)   
    check_fasta(genome)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        "simulate_variants",
        "conda" if not container else "conda apptainer",
        output_dir,
        2,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "simulate_variants.smk")

    conda_envs = ["simulations"]
    configs = {
        "workflow" : workflow,
        "prefix" : prefix,
        **({"random_seed" : random_seed} if random_seed else {}),
        "heterozygosity" : {
            "ratio" : heterozygosity,
            "only_vcf" : only_vcf,
        },
        "inversion" : {
            **({"vcf" : vcf} if vcf else {}),
            **({'count': count} if not vcf else {}),
            **({"min_size":  min_size} if not vcf else {}),
            **({"max_size" : max_size} if not vcf else {})
        },
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "inputs" : {
            "genome" : genome,
            **({"centromeres" : centromeres} if centromeres else {}),
            **({"genes" : genes} if genes else {}),
            **({"excluded_chromosomes" : exclude_chr} if exclude_chr else {})
        }
    }

    write_workflow_config(configs, output_dir)
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
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/simulate.inversion.summary")


@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Please Documentation: https://pdimens.github.io/harpy/workflows/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), help = 'VCF file of known copy number variants to simulate')
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
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Simulate/cnv", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prefix', type = str, default= "sim", show_default=True, help = "Naming prefix for output files")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
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
    Use `--only-vcf` alongside `--heterozygosity` to only generate the second VCF file and not simulate a second FASTA file.
 
    The two ratio parameters control different things and have special meanings when setting their value to either `9999` or `0`:
    
    | ratio | meaning | `9999` | `0` |
    |:----|:----|:---|:----|
    | `--dup-ratio` | tandem / dispersed | tand. only | disp. only |
    | `--gain-ratio` | copy gain / loss | gain only | loss only| 
    """
    workflow = "simulate_cnv"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow, True)
    ## checks and validations ##
    if not vcf and count == 0:
        print_error("missing option", "Provide either a `--count` of cnv to randomly simulate or a `--vcf` of known cnv to simulate.")
        sys.exit(1)
    if vcf and count > 0:
        print_error("conflicting arguments", "You can either simulate random CNVs (--count) or known CNVs (--vcf), but not both.")
        sys.exit(1) 
    check_fasta(genome)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        "simulate_variants",
        "conda" if not container else "conda apptainer",
        output_dir,
        2,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "simulate_variants.smk")

    conda_envs = ["simulations"]
    configs = {
        "workflow" : workflow,
        "prefix" : prefix,
        **({"random_seed" : random_seed} if random_seed else {}),
        "heterozygosity" : {
            "ratio" : heterozygosity,
            "only_vcf" : only_vcf,
        },
        "cnv" : {
            **({"vcf" : vcf} if vcf else {}),
            **({'count': count} if not vcf else {}),
            **({"min_size":  min_size} if not vcf else {}),
            **({"max_size" : max_size} if not vcf else {}),
            **({"duplication_ratio" : dup_ratio} if not vcf else {}),
            **({"max_copy" : max_copy} if not vcf else {}),
            **({"gain_ratio" : gain_ratio} if not vcf else {})
        },
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "inputs" : {
            "genome" : genome,
            **({"centromeres" : centromeres} if centromeres else {}),
            **({"genes" : genes} if genes else {}),
            **({"excluded_chromosomes" : exclude_chr} if exclude_chr else {})
        }
    }

    write_workflow_config(configs, output_dir)
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
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/simulate.cnv.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Please Documentation: https://pdimens.github.io/harpy/workflows/simulate/simulate-variants")
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), help = 'VCF file of known translocations to simulate')
@click.option('-n', '--count', type = click.IntRange(min = 0), default=0, show_default=False, help = "Number of random translocations to simluate")
@click.option('-c', '--centromeres', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of centromeres to avoid")
@click.option('-g', '--genes', type = InputFile("gff", gzip_ok = True), help = "GFF3 file of genes to avoid when simulating")
@click.option('-z', '--heterozygosity', type = click.FloatRange(0,1), default = 0, show_default=True, help = 'heterozygosity to simulate diploid variants')
@click.option('--only-vcf',  is_flag = True, default = False, help = 'If setting heterozygosity, only create the vcf rather than the fasta files')
@click.option('-e', '--exclude-chr', type = click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), help = "Text file of chromosomes to avoid")
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Simulate/translocation", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prefix', type = str, default= "sim", show_default=True, help = "Naming prefix for output files")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
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
    Use `--only-vcf` alongside `--heterozygosity` to only generate the second VCF file and not simulate a second FASTA file.
    """
    workflow = "simulate_translocation"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow, True)
    ## checks and validations ##
    if not vcf and count == 0:
        print_error("missing option", "Provide either a `--count` of cnv to randomly simulate or a `--vcf` of known cnv to simulate.")
        sys.exit(1)
    if vcf and count > 0:
        print_error("conflicting arguments", "You can either simulate random translocations (--count) or known translocations (--vcf), but not both.")
        sys.exit(1) 
    check_fasta(genome)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        "simulate_variants",
        "conda" if not container else "conda apptainer",
        output_dir,
        2,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "simulate_variants.smk")

    conda_envs = ["simulations"]
    configs = {
        "workflow" : workflow,
        "prefix" : prefix,
        **({"random_seed" : random_seed} if random_seed else {}),
        "heterozygosity" : {
            "ratio" : heterozygosity,
            "only_vcf" : only_vcf
        },
        "translocation" : {
            **({"vcf" : vcf} if vcf else {}),
            **({'count': count} if not vcf else {})
        },
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "inputs" : {
            "genome" : genome,
            **({"centromeres" : centromeres} if centromeres else {}),
            **({"genes" : genes} if genes else {}),
            **({"excluded_chromosomes" : exclude_chr} if exclude_chr else {})
        }
    }

    write_workflow_config(configs, output_dir)
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
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/simulate.translocation.summary")

