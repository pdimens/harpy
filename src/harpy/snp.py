from .helperfunctions import fetch_file, generate_conda_deps, getnames, createregions, print_onstart
from .helperfunctions import validate_bamfiles, validate_popfile, validate_vcfsamples
import rich_click as click
import subprocess
import sys
import os

@click.command(no_args_is_help = True)
@click.option('-g', '--genome', type=click.Path(exists=True), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-p', '--populations', type=click.Path(exists = True), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, metavar = "Integer", help = 'Ploidy of samples')
@click.option('-w', '--windowsize', default = 50000, show_default = True, type = int, metavar = "Integer", help = "Interval size for parallel variant calling")
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def mpileup(genome, threads, directory, populations, ploidy, windowsize, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call variants from using bcftools mpileup
    
    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    fetch_file("snp-mpileup.smk", "Variants/mpileup/workflow/")
    fetch_file("BcftoolsStats.Rmd", "Variants/mpileup/workflow/report/")

    samplenames = getnames(directory, '.bam')
    callcoords, linkedgenome = createregions(genome, windowsize, "mpileup")
    directory = directory.rstrip("/^")
    validate_bamfiles(directory, samplenames)
    command = (f'snakemake --rerun-incomplete --nolock --cores {threads} --directory .').split()
    command.append('--snakefile')
    command.append('Variants/mpileup/workflow/snp-mpileup.smk')
    command.append('--configfile')
    command.append('Variants/mpileup/workflow/config.yml')
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]

    call_SM = " ".join(command)

    with open("Variants/mpileup/workflow/config.yml", "w") as config:
        config.write(f"seq_directory: {directory}\n")
        config.write(f"samplenames: {samplenames}\n")
        popgroupings = ""
        if populations is not None:
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(directory, populations, samplenames, rows, quiet)
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

    if print_only:
        click.echo(call_SM)
    else:
        print_onstart(
            f"Initializing the [bold]harpy snp mpileup[/bold] workflow.\nInput Directory: {directory}\nSamples: {len(samplenames)}{popgroupings}"
        )
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)

@click.command(no_args_is_help = True)
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, metavar = "Integer", help = 'Ploidy of samples')
@click.option('-w', '--windowsize', default = 50000, show_default = True, type = int, metavar = "Integer", help = "Interval size for parallel variant calling")
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def freebayes(genome, threads, directory, populations, ploidy, windowsize, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call variants using freebayes
    
    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    fetch_file("snp-freebayes.smk", "Variants/freebayes/workflow/")
    fetch_file("BcftoolsStats.Rmd", "Variants/freebayes/workflow/report/")

    samplenames = getnames(directory, '.bam')
    callcoords, linkedgenome = createregions(genome, windowsize, "freebayes")
    directory = directory.rstrip("/^")
    validate_bamfiles(directory, samplenames)
    command = (f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake --cores {threads} --directory .').split()
    command.append('--snakefile')
    command.append('Variants/freebayes/workflow/snp-freebayes.smk')
    command.append('--configfile')
    command.append('Variants/freebayes/workflow/config.yml')
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]

    call_SM = " ".join(command)

    with open("Variants/mpileup/workflow/config.yml", "w") as config:
        config.write(f"seq_directory: {directory}\n")
        config.write(f"samplenames: {samplenames}\n")
        popgroupings = ""
        if populations is not None:
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(directory, populations, samplenames, rows, quiet)
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

    if print_only:
        click.echo(call_SM)
    else:
        print_onstart(
            f"Initializing the [bold]harpy snp freebayes[/bold] workflow.\nInput Directory: {directory}\nSamples: {len(samplenames)}{popgroupings}"
        )
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)