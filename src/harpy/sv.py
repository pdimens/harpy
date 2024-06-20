"""Harpy workflows to detect structural variants"""

import os
import sys
import subprocess
from pathlib import Path
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, fetch_report
from .fileparsers import parse_alignment_inputs
from .printfunctions import print_onstart
from .validations import validate_popfile, validate_popsamples, check_phase_vcf, validate_input_by_ext

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def sv():
    """
    Call large structural variants on alignments
 
    **Structural Variant Callers**
    - `naibr`: calls inversions, duplicates, deletions
    - `leviathan`: calls inversions, duplicates, deletions, misc breakends

    Provide an additional subcommand `leviathan` or `naibr` to get more information on using
    those variant callers. NAIBR tends to call variants better, but requires more user preprocessing.
    """

docstring = {
    "harpy sv leviathan": [
        {
            "name": "Parameters",
            "options": ["--genome", "--min-sv", "--min-barcodes", "--populations", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy sv naibr": [
        {
            "name": "Module Parameters",
            "options": ["--genome", "--vcf", "--min-sv", "--min-barcodes", "--molecule-distance", "--populations", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog= "read the docs for more information: https://pdimens.github.io/harpy/modules/sv/leviathan/")
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, help = 'Genome assembly for variant calling')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-n', '--min-sv', type = click.IntRange(min = 10, max_open = True), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-b', '--min-barcodes', show_default = True, default=2, type = click.IntRange(min = 1, max_open = True), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-o', '--output-dir', type = str, default = "SV/leviathan", show_default=True, help = 'Name of output directory')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def leviathan(inputs, output_dir, genome, min_sv, min_barcodes, threads, populations, extra_params, snakemake, skipreports, quiet, hpc, conda, print_only):
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
    sdm = "conda" if conda else "conda apptainer"
    vcaller = "leviathan" if populations is None else "leviathan-pop"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/sv-{vcaller}.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    if print_only:
        click.echo(command)
        sys.exit(0)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    if populations is not None:
        fetch_report(workflowdir, "LeviathanPop.Rmd")
    fetch_report(workflowdir, "Leviathan.Rmd")
    fetch_rule(workflowdir, f"sv-{vcaller}.smk")

    with open(f'{workflowdir}/config.yaml', "w", encoding="utf-8") as config:
        config.write("workflow: sv leviathan\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"min_barcodes: {min_barcodes}\n")
        config.write(f"min_sv: {min_sv}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        popgroupings = ""
        config.write(f"  genome: {Path(genome).resolve()}\n")
        if populations is not None:
            # check for delimeter and formatting
            validate_popfile(populations)
            # check that samplenames and populations line up
            validate_popsamples(bamlist, populations,quiet)
            config.write(f"  groupings: {Path(populations).resolve()}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write("  alignments:\n")
        for i in bamlist:
            config.write(f"    - {i}\n")

    modetext = "pool-by-group" if populations else "single-sample"
    print_onstart(
        f"Samples: {n}{popgroupings}\nOutput Directory: {output_dir}/\nMode: {modetext}",
        "sv leviathan"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

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
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def naibr(inputs, output_dir, genome, vcf, min_sv, min_barcodes, threads, populations, molecule_distance, extra_params, snakemake, skipreports, quiet, hpc, conda, print_only):
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
    sdm = "conda" if conda else "conda apptainer"
    vcaller = "naibr" if populations is None else "naibr-pop"
    vcaller += "-phase" if vcf is not None else ""
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/sv-{vcaller}.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake
    if print_only:
        click.echo(command)
        sys.exit(0)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_input_by_ext(genome, "--genome", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    if populations is not None:
        fetch_report(workflowdir, "NaibrPop.Rmd")
    fetch_report(workflowdir, "Naibr.Rmd")
    if vcf is not None:
        check_phase_vcf(vcf)
        #vcaller += "-phase"
    fetch_rule(workflowdir, f"sv-{vcaller}.smk")

    with open(f'{workflowdir}/config.yaml', "w", encoding="utf-8") as config:
        config.write("workflow: sv naibr\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"min_barcodes: {min_barcodes}\n")
        config.write(f"min_sv: {min_sv}\n")
        config.write(f"molecule_distance: {molecule_distance}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        popgroupings = ""
        config.write("inputs:\n")
        if vcf is not None:
            config.write(f"  vcf: {Path(vcf).resolve()}\n")
        if genome is not None:
            config.write(f"  genome: {Path(genome).resolve()}\n")
        if populations is not None:
            # check for delimeter and formatting
            validate_popfile(populations)
            # check that samplenames and populations line up
            validate_popsamples(bamlist, populations, quiet)
            config.write(f"  groupings: {Path(populations).resolve()}\n")
            popgroupings += f"\nPopulations: {populations}"
        config.write("  alignments:\n")
        for i in bamlist:
            config.write(f"    - {i}\n")

    if populations:
        modetext = "pool-by-group"
    else:
        modetext = "single-sample"
    if vcf:
        modetext += " + will be phased"
    else:
        modetext += " + already phased"
    print_onstart(
        f"Samples: {n}{popgroupings}\nOutput Directory: {output_dir}/\nMode: {modetext}",
        "sv naibr"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)


sv.add_command(leviathan)
sv.add_command(naibr)
