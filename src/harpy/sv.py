from .helperfunctions import fetch_rule, fetch_report, generate_conda_deps
from .fileparsers import getnames, parse_alignment_inputs
from .printfunctions import print_onstart
from .validations import validate_popfile, validate_vcfsamples, check_phase_vcf, validate_input_by_ext
import rich_click as click
import subprocess
import sys
import os

@click.command(no_args_is_help = True, epilog= "read the docs for more information: https://pdimens.github.io/harpy/modules/sv/leviathan/")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), metavar = "File Path", help = "Tab-delimited file of sample\<tab\>population")
@click.option('-n', '--min-sv', type = click.IntRange(min = 10, max_open = True), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-b', '--min-barcodes', show_default = True, default=2, type = click.IntRange(min = 1, max_open = True), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-o', '--output-dir', type = str, default = "SV/leviathan", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def leviathan(input, output_dir, genome, min_sv, min_barcodes, threads, populations, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Call structural variants using LEVIATHAN
    
    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    Optionally specify `--populations` for population-pooled variant calling. 
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    vcaller = "leviathan" if populations is None else "leviathan-pop"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method conda apptainer --use-apptainer --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
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
    validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    if populations is not None:
        fetch_report(workflowdir, "LeviathanPop.Rmd")
    fetch_report(workflowdir, "Leviathan.Rmd")
    fetch_rule(workflowdir, f"sv-{vcaller}.smk")

    with open(f'{workflowdir}/config.yml', "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
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
        config.write(f"min_barcodes: {min_barcodes}\n")
        config.write(f"min_sv: {min_sv}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Samples: {len(samplenames)}{popgroupings}\nOutput Directory: {output_dir}/",
        "sv leviathan"
    )
    return command

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/sv/naibr/")
@click.option('-g', '--genome', required = True, type=click.Path(exists=True, dir_okay=False), help = 'Genome assembly')
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False),  help = 'Path to phased bcf/vcf file')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-n', '--min-sv', type = click.IntRange(min = 10, max_open = True), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-b', '--min-barcodes', show_default = True, default=2, type = click.IntRange(min = 1, max_open = True), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Base-pair distance delineating separate molecules')
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-o', '--output-dir', type = str, default = "SV/naibr", show_default=True, help = 'Name of output directory')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate any HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def naibr(input, output_dir, genome, vcf, min_sv, min_barcodes, threads, populations, molecule_distance, extra_params, snakemake, skipreports, quiet, print_only):
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
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    vcaller = "naibr" if populations is None else "naibr-pop"
    vcaller += "-phase" if vcf is not None else ""
    
    command = (f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method conda apptainer --use-apptainer --conda-prefix ./.snakemake/conda --cores {threads} --directory .').split()
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
    validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    if populations is not None:
        fetch_report(workflowdir, "NaibrPop.Rmd")
    fetch_report(workflowdir, "Naibr.Rmd")
    if vcf is not None:
        check_phase_vcf(vcf)
        #vcaller += "-phase"
    fetch_rule(workflowdir, f"sv-{vcaller}.smk")

    with open(f'{workflowdir}/config.yml', "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
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
        config.write(f"min_barcodes: {min_barcodes}\n")
        config.write(f"min_sv: {min_sv}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Samples: {len(samplenames)}{popgroupings}\nOutput Directory: {output_dir}/",
        "sv naibr"
    )
    return command
