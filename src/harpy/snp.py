from .helperfunctions import fetch_rule, fetch_report, fetch_script, generate_conda_deps, createregions
from .fileparsers import getnames, parse_alignment_inputs
from .printfunctions import print_onstart
from .validations import validate_bamfiles, validate_popfile, validate_vcfsamples, validate_input_by_ext
import rich_click as click
import sys
import os

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/snp")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), metavar = "File Path", help = "Tab-delimited file of sample\<tab\>population")
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, metavar = "Integer", help = 'Ploidy of samples')
@click.option('-w', '--windowsize', default = 50000, show_default = True, type = int, metavar = "Integer", help = "Interval size for parallel variant calling")
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "SNP/mpileup", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def mpileup(input, output_dir, genome, threads, populations, ploidy, windowsize, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call variants from using bcftools mpileup
    
    Provide the input alignment (`.bam`) files and/or directories
    at the end of the command as individual files/folders, using shell wildcards
    (e.g. `data/scarab*.bam`), or both.
    
    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = (f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .').split()
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
    fetch_rule(workflowdir, "snp-mpileup.smk")
    fetch_report(workflowdir, "BcftoolsStats.Rmd")
    callcoords, linkedgenome = createregions(genome, windowsize, "mpileup")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplenames: {samplenames}\n")
        popgroupings = ""
        if populations is not None:
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(f"{workflowdir}/input", populations, samplenames, rows, quiet)
            config.write(f"groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write(f"genomefile: {linkedgenome}\n")
        config.write(f"ploidy: {ploidy}\n")
        config.write(f"windowsize: {windowsize}\n")
        config.write(f"intervals: {callcoords}\n")
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
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), metavar = "File Path", help = "Tab-delimited file of sample\<tab\>population")
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, metavar = "Integer", help = 'Ploidy of samples')
@click.option('-w', '--windowsize', default = 50000, show_default = True, type = int, metavar = "Integer", help = "Interval size for parallel variant calling")
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "SNP/freebayes", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def freebayes(input, output_dir, genome, threads, populations, ploidy, windowsize, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call variants using freebayes
    
    Provide the input alignment (`.bam`) files and/or directories
    at the end of the command as individual files/folders, using shell wildcards
    (e.g. `data/jellyfish*.bam`), or both.
    
    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = (f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .').split()
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
    fetch_rule(workflowdir, "snp-freebayes.smk")
    fetch_report(workflowdir, "BcftoolsStats.Rmd")
    callcoords, linkedgenome = createregions(genome, windowsize, "freebayes")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplenames: {samplenames}\n")
        popgroupings = ""
        if populations is not None:
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(f"{workflowdir}/input", populations, samplenames, rows, quiet)
            config.write(f"groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write(f"ploidy: {ploidy}\n")
        config.write(f"windowsize: {windowsize}\n")
        config.write(f"intervals: {callcoords}\n")
        config.write(f"genomefile: {linkedgenome}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Samples: {len(samplenames)}{popgroupings}\nOutput Directory: {output_dir}/",
        "snp freebayes"
    )
    return command