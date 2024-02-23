from .helperfunctions import fetch_file, generate_conda_deps, getnames, print_onstart
from .helperfunctions import validate_popfile, validate_vcfsamples, check_phase_vcf, parse_alignment_inputs
import rich_click as click
import subprocess
import sys
import os

@click.command(no_args_is_help = True, epilog= "read the docs for more information: https://pdimens.github.io/harpy/modules/sv/leviathan/")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def leviathan(input, genome, threads, populations, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call structural variants using LEVIATHAN
    
    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    Optionally specify `--populations` for population-pooled variant calling. 
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    vcaller = "leviathan" if populations is None else "leviathan-pop"
    workflowdir = f"Variants/{vcaller}/workflow"

    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/sv-{vcaller}.smk')
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
    _ = parse_alignment_inputs(input, f"{workflowdir}/input")
    samplenames = getnames(f"{workflowdir}/input", '.bam')
    if populations is not None:
        fetch_file("LeviathanPop.Rmd", f"{workflowdir}/report/")
    fetch_file("Leviathan.Rmd", f"{workflowdir}/report/")
    fetch_file(f"sv-{vcaller}.smk", f"{workflowdir}")

    with open(f'{workflowdir}/config.yml', "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"samplenames: {samplenames}\n")
        popgroupings = ""
        if populations is not None:
            # check for delimeter and formatting
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(f"{workflowdir}/input", populations, samplenames, rows, quiet)
            config.write(f"groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write(f"genomefile: {genome}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        f"Samples: {len(samplenames)}{popgroupings}\nOutput Directory: Variants/{vcaller}/",
        "sv leviathan"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/sv/naibr/")
@click.option('-g', '--genome', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Genome assembly')
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Path to phased bcf/vcf file')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance delineating separate molecules')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def naibr(input, genome, vcf, threads, populations, molecule_distance, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call structural variants using NAIBR
    
    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    NAIBR requires **phased** bam files as input. This appears as the `HP` or `PS` tags
    in alignment records. If your bam files do not have either of these phasing tags
    (e.g. BWA/EMA do not phase alignments), then provide a **phased** `--vcf` file such
     as that created by `harpy phase` and Harpy will use [whatshap haplotag](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)
    to phase your input bam files prior to calling variants with NAIBR.

    Optionally specify `--populations` for population-pooled variant calling. 
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    vcaller = "naibr" if populations is None else "naibr-pop"
    workflowdir = f"Variants/{vcaller}/workflow"
    vcaller += "-phase" if vcf is not None else ""
    
    command = (f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .').split()
    command.append('--snakefile')
    command.append(f"{workflowdir}/sv-{vcaller}.smk")
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

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    _ = parse_alignment_inputs(input, f"{workflowdir}/input")
    samplenames = getnames(f"{workflowdir}/input", '.bam')
    if populations is not None:
        fetch_file("NaibrPop.Rmd", f"{workflowdir}/report/")
    fetch_file("Naibr.Rmd", f"{workflowdir}/report/")
    if vcf is not None:
        check_phase_vcf(vcf)
        #vcaller += "-phase"
    fetch_file(f"sv-{vcaller}.smk", f"{workflowdir}/")

    with open(f'{workflowdir}/config.yml', "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"samplenames: {samplenames}\n")
        popgroupings = ""
        if populations is not None:
            # check for delimeter and formatting
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(f"{workflowdir}/input", populations, samplenames, rows, quiet)
            config.write(f"groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write(f"molecule_distance: {molecule_distance}\n")
        if vcf is not None:
            config.write(f"vcf: {vcf}\n")
        if genome is not None:
            config.write(f"genomefile: {genome}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        f"Samples: {len(samplenames)}{popgroupings}\nOutput Directory: Variants/{vcaller}/",
        "sv naibr"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)