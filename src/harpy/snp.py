from .helperfunctions import fetch_rule, fetch_report, fetch_script, generate_conda_deps, createregions
from .fileparsers import getnames, parse_alignment_inputs
from .printfunctions import print_onstart, print_error
from .validations import validate_bamfiles, validate_popfile, validate_vcfsamples, validate_input_by_ext
import rich_click as click
import subprocess
import sys
import os

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/snp")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, help = 'Ploidy of samples')
@click.option('-w', '--windowsize', type = click.IntRange(min = 10), help = "Interval size for parallel variant calling (default: 50000)")
@click.option('-r', '--regions', type=click.Path(exists = True, dir_okay=False), help = "BED or tab-delimited text file of regions to call variants")
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "SNP/mpileup", show_default=True, help = 'Name of output directory')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate any HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def mpileup(input, output_dir, regions, genome, threads, populations, ploidy, windowsize, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call variants from using bcftools mpileup
    
    Provide the input alignment (`.bam`) files and/or directories
    at the end of the command as individual files/folders, using shell wildcards
    (e.g. `data/scarab*.bam`), or both.
    
    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`. 
    """
    if regions and windowsize:
        print_error("The options [yellow bold]--regions[/yellow bold] and [yellow bold]--windowsize[/yellow bold] cannot be used together")
        exit(1)
    
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
    if windowsize>0:
        window_size = 50000 if windowsize==0 else windowsize
        linkedgenome = createregions(genome, window_size, "mpileup", f"{workflowdir}/positions.bed")
    else:
        bn = os.path.basename(genome)
        ftype = subprocess.run(["file", genome], stdout=subprocess.PIPE).stdout.decode('utf-8')
        if "Blocked GNU Zip" in ftype:
            # is bgzipped, just link it
            subprocess.run(f"ln -sr {genome} Genome/{bn}".split())
        elif "gzip compressed data" in ftype:
            # is regular gzipped, needs to be bgzipped
            subprocess.run(f"zcat {genome} | bgzip -c > Genome/{bn}".split())
        else:
            # not compressed, just link
            subprocess.run(f"ln -sr {genome} Genome/{bn}".split())
        # TODO MAKE SURE THIS FOLLOWS THE f"{contig}:{startpos}-{endpos}\n" FORMAT
        os.system(f"cp -f {regions} {workflowdir}/positions.bed")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplenames: {samplenames}\n")
        if regions:
            config.write(f"regions: {regions}\n")
        else:
            config.write(f"windowsize: {window_size}\n")
            config.write(f"intervals: {callcoords}\n")
        popgroupings = ""
        if populations is not None:
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(f"{workflowdir}/input", populations, samplenames, rows, quiet)
            config.write(f"groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write(f"genomefile: {linkedgenome}\n")
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
@click.option('-w', '--windowsize', default = 50000, show_default = True, type = int, help = "Interval size for parallel variant calling")
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "SNP/freebayes", show_default=True, help = 'Name of output directory')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate any HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
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