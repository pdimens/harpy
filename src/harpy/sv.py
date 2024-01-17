from .helperfunctions import fetch_file, generate_conda_deps, getnames, check_phase_vcf
from .helperfunctions import validate_popfile, validate_vcfsamples
import rich_click as click
import subprocess
import sys
import os

@click.command(no_args_is_help = True)
@click.option('-g', '--genome', type=click.Path(exists=True), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-p', '--populations', type=click.Path(exists = True), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def leviathan(genome, threads, directory, populations, extra_params, snakemake, quiet, print_only):
    """
    Call structural variants using LEVIATHAN
    
    Optionally specify `--populations` for population-pooled variant calling. 
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    samplenames = getnames(directory, '.bam')
    vcaller = "leviathan"
    if populations is not None:
        vcaller += "-pop"
        fetch_file("LeviathanPop.Rmd", f"Variants/{vcaller}/workflow/report/")
    fetch_file("Leviathan.Rmd", f"Variants/{vcaller}/workflow/report/")
    fetch_file(f"variants-{vcaller}.smk", f"Variants/{vcaller}/workflow/")
    directory = directory.rstrip("/^")
    command = (f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake --cores {threads} --directory . --snakefile Variants/{vcaller}/workflow/variants-{vcaller}.smk').split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"seq_directory={directory}")
    command.append(f"samplenames={samplenames}")
    if populations is not None:
        # check for delimeter and formatting
        rows = validate_popfile(populations)
        # check that samplenames and populations line up
        validate_vcfsamples(directory, populations, samplenames, rows, quiet)
        command.append(f"groupings={populations}")
    command.append(f"genomefile={genome}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    if print_only:
        click.echo(" ".join(command))
    else:
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)

@click.command(no_args_is_help = True)
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-g', '--genome', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Genome assembly')
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Path to phased bcf/vcf file')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance delineating separate molecules')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def naibr(genome, vcf, threads, directory, populations, molecule_distance, extra_params, snakemake, quiet, print_only):
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
    samplenames = getnames(directory, '.bam')
    vcaller = "naibr"
    outdir = "naibr"
    if populations is not None:
        vcaller += "-pop"
        outdir += "-pop"
        fetch_file("NaibrPop.Rmd", f"Variants/{vcaller}/workflow/report/")
    fetch_file("Naibr.Rmd", f"Variants/{vcaller}/workflow/report/")
    if vcf is not None:
        check_phase_vcf(vcf)
        vcaller += "-phase"
    fetch_file(f"variants-{vcaller}.smk", f"Variants/{vcaller}/workflow/")
    command = (f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake --cores {threads} --directory . --snakefile Variants/{outdir}/workflow/variants-{vcaller}.smk').split()
    
    directory = directory.rstrip("/^")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"seq_directory={directory}")
    command.append(f"samplenames={samplenames}")
    if populations is not None:
        # check for delimeter and formatting
        rows = validate_popfile(populations)
        # check that samplenames and populations line up
        validate_vcfsamples(directory, populations, samplenames, rows, quiet)
        command.append(f"groupings={populations}")
    command.append(f"molecule_distance={molecule_distance}")
    if vcf is not None:
        command.append(f"vcf={vcf}")
    if genome is not None:
        command.append(f"genomefile={genome}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    if print_only:
        click.echo(" ".join(command))
    else:
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)