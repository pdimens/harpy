from .helperfunctions import fetch_file, generate_conda_deps, getnames, print_onstart
from .helperfunctions import validate_popfile, validate_vcfsamples, check_phase_vcf
import rich_click as click
import subprocess
import sys
import os

@click.command(no_args_is_help = True, epilog= "read the docs for more information: https://pdimens.github.io/harpy/modules/sv/leviathan/")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def leviathan(genome, threads, directory, populations, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call structural variants using LEVIATHAN
    
    Optionally specify `--populations` for population-pooled variant calling. 
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    directory = directory.rstrip("/^")
    samplenames = getnames(directory, '.bam')
    vcaller = "leviathan"
    if populations is not None:
        vcaller += "-pop"
        fetch_file("LeviathanPop.Rmd", f"Variants/{vcaller}/workflow/report/")
    fetch_file("Leviathan.Rmd", f"Variants/{vcaller}/workflow/report/")
    fetch_file(f"sv-{vcaller}.smk", f"Variants/{vcaller}/workflow/")
    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'Variants/{vcaller}/workflow/sv-{vcaller}.smk')
    command.append('--configfile')
    command.append(f'Variants/{vcaller}/workflow/config.yml')
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]

    call_SM = " ".join(command)

    with open(f'Variants/{vcaller}/workflow/config.yml', "w") as config:
        config.write(f"seq_directory: {directory}\n")
        config.write(f"samplenames: {samplenames}\n")
        popgroupings = ""
        if populations is not None:
            # check for delimeter and formatting
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(directory, populations, samplenames, rows, quiet)
            config.write(f"groupings: {populations}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write(f"genomefile: {genome}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    if print_only:
        click.echo(call_SM)
    else:
        print_onstart(
            f"Input Directory: {directory}\nSamples: {len(samplenames)}{popgroupings}",
            "sv leviathan"
        )
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/sv/naibr/")
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with BAM alignments')
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
def naibr(genome, vcf, threads, directory, populations, molecule_distance, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call structural variants using NAIBR
    
    NAIBR requires **phased** bam files as input. This appears as the `HP` or `PS` tags
    in alignment records. If your bam files do not have either of these phasing tags
    (e.g. BWA/EMA do not phase alignments), then provide a **phased** `--vcf` file such
     as that created by `harpy phase` and Harpy will use [whatshap haplotag](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)
    to phase your input bam files prior to calling variants with NAIBR.

    Optionally specify `--populations` for population-pooled variant calling. 
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    directory = directory.rstrip("/^")
    samplenames = getnames(directory, '.bam')
    vcaller = "naibr"
    outdir = "naibr"
    if populations is not None:
        vcaller += "-pop"
        outdir += "-pop"
        fetch_file("NaibrPop.Rmd", f"Variants/{outdir}/workflow/report/")
    fetch_file("Naibr.Rmd", f"Variants/{outdir}/workflow/report/")
    if vcf is not None:
        check_phase_vcf(vcf)
        vcaller += "-phase"
    fetch_file(f"sv-{vcaller}.smk", f"Variants/{outdir}/workflow/")
    command = (f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake --cores {threads} --directory .').split()
    command.append('--snakefile')
    command.append(f"Variants/{outdir}/workflow/sv-{vcaller}.smk")
    command.append("--configfile")
    command.append(f"Variants/{outdir}/workflow/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]

    call_SM = " ".join(command)

    with open(f'Variants/{outdir}/workflow/config.yml', "w") as config:
        config.write(f"seq_directory: {directory}\n")
        config.write(f"samplenames: {samplenames}\n")
        popgroupings = ""
        if populations is not None:
            # check for delimeter and formatting
            rows = validate_popfile(populations)
            # check that samplenames and populations line up
            validate_vcfsamples(directory, populations, samplenames, rows, quiet)
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

    if print_only:
        click.echo(call_SM)
    else:
        print_onstart(
            f"Input Directory: {directory}\nSamples: {len(samplenames)}{popgroupings}",
            "sv naibr"
        )
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)