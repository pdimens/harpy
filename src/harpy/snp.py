"""Harpy workflows to call SNP variants"""

import os
import sys
import subprocess
from pathlib import Path
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, fetch_report
from .fileparsers import parse_alignment_inputs
from .printfunctions import print_onstart
from .validations import validate_bam_RG, validate_popfile, validate_popsamples, validate_input_by_ext, validate_regions

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def snp():
    """
    Call SNPs and small indels on alignments
    
    **Variant Callers**
    - `mpileup`: call variants using bcftools mpileup
    - `freebayes`: call variants using freebayes

    Provide an additional subcommand `mpileup` or `freebayes` to get more information on using
    those variant callers. They are both robust variant callers and neither is recommended over the other.
    """

docstring = {
    "harpy snp mpileup": [
        {
            "name": "Parameters",
            "options": ["--genome", "--populations", "--ploidy", "--regions", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy snp freebayes": [
        {
            "name": "Parameters",
            "options": ["--genome", "--populations", "--ploidy", "--regions", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--hpc",  "--conda", "--snakemake", "--quiet", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/snp")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, help = 'Ploidy of samples')
@click.option('-r', '--regions', type=str, default=50000, show_default=True, help = "Regions where to call variants")
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-o', '--output-dir', type = str, default = "SNP/mpileup", show_default=True, help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--config-only',  is_flag = True, hidden = True, default = False, help = 'Create the config.yaml file and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def mpileup(inputs, output_dir, regions, genome, threads, populations, ploidy, extra_params, snakemake, skipreports, quiet, hpc, conda, config_only):
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
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/snp-mpileup.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_bam_RG(bamlist)
    validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])

    # setup regions checks
    regtype = validate_regions(regions, genome)
    if regtype == "windows":
        region = Path(f"{workflowdir}/positions.bed").resolve()
        os.system(f"makeWindows.py -m 1 -i {genome} -o {region} -w {regions}")
    elif regtype == "region":
        region = regions
    else:
        region = Path(f"{workflowdir}/positions.bed").resolve()
        os.system(f"cp -f {regions} {region}")

    fetch_rule(workflowdir, "snp-mpileup.smk")
    fetch_report(workflowdir, "BcftoolsStats.Rmd")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: snp mpileup\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"ploidy: {ploidy}\n")
        config.write(f"regiontype: {regtype}\n")
        if regtype == "windows":
            config.write("windowsize: {regions}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome: {Path(genome).resolve()}\n")
        config.write(f"  regions: {region}\n")
        popgroupings = ""
        if populations is not None:
            validate_popfile(populations)
            # check that samplenames and populations line up
            validate_popsamples(bamlist, populations, quiet)
            config.write(f"  groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write("  alignments:\n")
        for i in bamlist:
            config.write(f"    - {i}\n")
    if config_only:
        sys.exit(0)

    print_onstart(
        f"Samples: {n}{popgroupings}\nOutput Directory: {output_dir}/",
        "snp mpileup"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/snp")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, help = 'Ploidy of samples')
@click.option('-r', '--regions', type=str, default=50000, show_default=True, help = "Regions where to call variants")
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-o', '--output-dir', type = str, default = "SNP/freebayes", show_default=True, help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--config-only',  is_flag = True, hidden = True, default = False, help = 'Create the config.yaml file and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def freebayes(inputs, output_dir, genome, threads, populations, ploidy, regions, extra_params, snakemake, skipreports, quiet, hpc, conda, config_only):
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
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/snp-freebayes.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_bam_RG(bamlist)
    validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])

    # setup regions checks
    regtype = validate_regions(regions, genome)
    if regtype == "windows":
        region = Path(f"{workflowdir}/positions.bed").resolve()
        os.system(f"makeWindows.py -m 0 -i {genome} -o {region} -w {regions}")
    elif regtype == "region":
        region = regions
    else:
        region = Path(f"{workflowdir}/positions.bed").resolve()
        os.system(f"cp -f {regions} {region}")

    fetch_rule(workflowdir, "snp-freebayes.smk")
    fetch_report(workflowdir, "BcftoolsStats.Rmd")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: snp freebayes\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"ploidy: {ploidy}\n")
        config.write(f"regiontype: {regtype}\n")
        if regtype == "windows":
            config.write("windowsize: {regions}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  genome: {Path(genome).resolve()}\n")
        config.write(f"  regions: {region}\n")
        popgroupings = ""
        if populations is not None:
            validate_popfile(populations)
            # check that samplenames and populations line up
            validate_popsamples(bamlist, populations, quiet)
            config.write(f"  groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write("  alignments:\n")
        for i in bamlist:
            config.write(f"    - {i}\n")
    if config_only:
        sys.exit(0)

    print_onstart(
        f"Samples: {n}{popgroupings}\nOutput Directory: {output_dir}/",
        "snp freebayes"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

snp.add_command(mpileup)
snp.add_command(freebayes)