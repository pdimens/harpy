from .helperfunctions import fetch_rule, fetch_report, fetch_script, generate_conda_deps
from .fileparsers import getnames, parse_alignment_inputs
from .printfunctions import print_onstart, print_error
from .validations import validate_bamfiles, validate_popfile, validate_vcfsamples, validate_input_by_ext, validate_regions
import rich_click as click
import subprocess
import sys
import os

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/snp")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, help = 'Ploidy of samples')
@click.option('-r', '--regions', type=str, default=50000, show_default=True, help = "Regions where to call variants")
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "SNP/mpileup", show_default=True, help = 'Name of output directory')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate any HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def mpileup(input, output_dir, regions, genome, threads, populations, ploidy, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call variants from using bcftools mpileup
    
    Provide the input alignment (`.bam`) files and/or directories
    at the end of the command as individual files/folders, using shell wildcards
    (e.g. `data/scarab*.bam`), or both.
    
    The `--regions` option specifies what genomic regions to call variants
    with. If a BED or tab delimited file is provided, variant calling will be parallelized
    over those regions. If a single region is provided in the format `chrom:start-end`, only
    that region will be called. If an integer is provided (default), then Harpy will
    call variants in parallel for intervals of that size across the entire genome.

    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`. 
    """   
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = (f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .').split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/snp-mpileup.smk')
    command.append('--configfile')
    command.append(f'{workflowdir}/config.yml')
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
    samplenames = getnames(f"{workflowdir}/input", '.bam')
    validate_bamfiles(f"{workflowdir}/input", samplenames)
    validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])

    # setup regions checks
    regtype = validate_regions(regions, genome)
    if regtype == "windows":
        region = f"{workflowdir}/positions.bed"
        os.system(f"makeWindows.py -m 1 -i {genome} -o {region} -w {regions}")
    elif regtype == "region":
        region = regions
    else:
        region = f"{workflowdir}/positions.bed"
        os.system(f"cp -f {regions} {region}")

    fetch_rule(workflowdir, "snp-mpileup.smk")
    fetch_report(workflowdir, "BcftoolsStats.Rmd")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplenames: {samplenames}\n")
        config.write(f"regiontype: {regtype}\n")
        if regtype == "windows":
            config.write("windowsize: {regions}\n")
        config.write(f"regions: {region}\n")
        popgroupings = ""
        if populations is not None:
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(f"{workflowdir}/input", populations, samplenames, rows, quiet)
            config.write(f"groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write(f"genomefile: {genome}\n")
        config.write(f"ploidy: {ploidy}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Samples: {len(samplenames)}{popgroupings}\nOutput Directory: {output_dir}/",
        "snp mpileup"
    )
    return command

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/snp")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, help = 'Ploidy of samples')
@click.option('-r', '--regions', type=str, default=50000, show_default=True, help = "Regions where to call variants")
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "SNP/freebayes", show_default=True, help = 'Name of output directory')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate any HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def freebayes(input, output_dir, genome, threads, populations, ploidy, regions, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call variants using freebayes
    
    Provide the input alignment (`.bam`) files and/or directories
    at the end of the command as individual files/folders, using shell wildcards
    (e.g. `data/jellyfish*.bam`), or both.
    
    The `--regions` option specifies what genomic regions to call variants
    with. If a BED or tab delimited file is provided, variant calling will be parallelized
    over those regions. If a single region is provided in the format `chrom:start-end`, only
    that region will be called. If an integer is provided (default), then Harpy will
    call variants in parallel for intervals of that size across the entire genome.

    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = (f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .').split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/snp-freebayes.smk')
    command.append('--configfile')
    command.append(f'{workflowdir}/config.yml')
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
    samplenames = getnames(f"{workflowdir}/input", '.bam')
    validate_bamfiles(f"{workflowdir}/input", samplenames)
    validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])

    # setup regions checks
    regtype = validate_regions(regions, genome)
    if regtype == "windows":
        region = f"{workflowdir}/positions.bed"
        os.system(f"makeWindows.py -m 0 -i {genome} -o {region} -w {regions}")
    elif regtype == "region":
        region = regions
    else:
        region = f"{workflowdir}/positions.bed"
        os.system(f"cp -f {regions} {region}")

    fetch_rule(workflowdir, "snp-freebayes.smk")
    fetch_report(workflowdir, "BcftoolsStats.Rmd")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplenames: {samplenames}\n")
        config.write(f"regiontype: {regtype}\n")
        if regtype == "windows":
            config.write("windowsize: {regions}\n")
        config.write(f"regions: {region}\n")
        popgroupings = ""
        if populations is not None:
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(f"{workflowdir}/input", populations, samplenames, rows, quiet)
            config.write(f"groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write(f"ploidy: {ploidy}\n")
        config.write(f"genomefile: {genome}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Samples: {len(samplenames)}{popgroupings}\nOutput Directory: {output_dir}/",
        "snp freebayes"
    )
    return command