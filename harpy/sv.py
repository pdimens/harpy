"""Harpy workflows to detect structural variants"""

import os
import sys
import yaml
from pathlib import Path
from rich import box
from rich.table import Table
import rich_click as click
from ._cli_types_generic import ContigList, HPCProfile, InputFile, SnakemakeParams
from ._cli_types_params import LeviathanParams, NaibrParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, fetch_report, snakemake_log
from ._parsers import parse_alignment_inputs
from ._validations import check_fasta, check_phase_vcf
from ._validations import validate_popfile, validate_popsamples, fasta_contig_match

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
            "options": ["--extra-params", "--genome", "--iterations", "--min-barcodes", "--min-sv", "--populations"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--contigs", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ],
    "harpy sv naibr": [
        {
            "name": "Module Parameters",
            "options": ["--extra-params", "--genome", "--min-barcodes", "--min-quality", "--min-sv", "--molecule-distance", "--populations", "--vcf"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--contigs", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/sv/leviathan/")
@click.option('-x', '--extra-params', type = str, help = 'Additional leviathan parameters, in quotes')
@click.option('-g', '--genome', type=InputFile("fasta", gzip_ok = True), required = True, help = 'Genome assembly for variant calling')
@click.option('-i', '--iterations', show_default = True, default=50, type = click.IntRange(min = 10, max_open = True), help = 'Number of iterations to perform through index (reduces memory)')
@click.option('-s', '--min-sv', type = click.IntRange(min = 10, max_open = True), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-b', '--min-barcodes', show_default = True, default=2, type = click.IntRange(min = 1, max_open = True), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "SV/leviathan", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of a container')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def leviathan(inputs, output_dir, genome, min_sv, min_barcodes, iterations, threads, populations, extra_params, snakemake, skip_reports, quiet, hpc, conda, contigs, setup_only):
    """
    Call structural variants using LEVIATHAN
    
    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    Optionally specify `--populations` for population-pooled variant calling. 
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if conda else "conda apptainer"
    vcaller = "leviathan" if populations is None else "leviathan_pop"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/sv_{vcaller}.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    check_fasta(genome)
    if contigs:
        fasta_contig_match(contigs, genome)
    if populations:
        validate_popfile(populations)
        validate_popsamples(bamlist, populations,quiet)
        fetch_report(workflowdir, "leviathan_pop.Rmd")
    fetch_report(workflowdir, "leviathan.Rmd")
    fetch_rule(workflowdir, f"sv_{vcaller}.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "sv_leviathan")
    conda_envs = ["align", "r", "variants"]
    configs = {
        "workflow" : "sv leviathan",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "min_barcodes" : min_barcodes,
        "min_sv" : min_sv,
        "iterations" : iterations,
        **({'extra': extra_params} if extra_params else {}),
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "genome" : Path(genome).resolve().as_posix(),
            **({'groupings': Path(populations).resolve().as_posix()} if populations else {}),
            "alignments" : [i.as_posix() for i in bamlist]
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column(header="value", justify="left")
    start_text.add_row("Samples:", f"{n}")
    start_text.add_row("Genome:", genome)
    start_text.add_row("Sample Pooling:", populations if populations else "no")
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, "sv_leviathan", start_text, output_dir, sm_log, quiet, "workflow/sv.leviathan.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/sv/naibr/")
@click.option('-x', '--extra-params', type = str, help = 'Additional naibr parameters, in quotes')
@click.option('-g', '--genome', type=InputFile("fasta", gzip_ok = True), required = True, help = 'Genome assembly for calling variants')
@click.option('-b', '--min-barcodes', show_default = True, default=2, type = click.IntRange(min = 1, max_open = True), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-q', '--min-quality', show_default = True, default=30, type = click.IntRange(min = 0, max = 40), help = 'Minimum mapping quality of reads to use')
@click.option('-s', '--min-sv', type = click.IntRange(min = 10, max_open = True), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Base-pair distance delineating separate molecules')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "SV/naibr", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True),  help = 'Path to phased bcf/vcf file')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of a container')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def naibr(inputs, output_dir, genome, vcf, min_sv, min_barcodes, min_quality, threads, populations, molecule_distance, extra_params, snakemake, skip_reports, quiet, hpc, conda, contigs, setup_only):
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
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if conda else "conda apptainer"
    vcaller = "naibr" if populations is None else "naibr_pop"
    vcaller += "_phase" if vcf else ""
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/sv_{vcaller}.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    check_fasta(genome)
    if contigs:
        fasta_contig_match(contigs, genome)
    if populations:
        validate_popfile(populations)
        validate_popsamples(bamlist, populations, quiet)
        fetch_report(workflowdir, "naibr_pop.Rmd")
    fetch_report(workflowdir, "naibr.Rmd")
    if vcf:
        check_phase_vcf(vcf)
    fetch_rule(workflowdir, f"sv_{vcaller}.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "sv_naibr")
    conda_envs = ["phase", "r", "variants"]
    configs = {
        "workflow" : "sv naibr",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "min_barcodes" : min_barcodes,
        "min_quality" : min_quality,
        "min_sv" : min_sv,
        "molecule_distance" : molecule_distance,
        **({'extra': extra_params} if extra_params else {}),
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            **({'genome': Path(genome).resolve().as_posix()} if genome else {}),
            **({'vcf': Path(vcf).resolve().as_posix()} if vcf else {}),
            **({'groupings': Path(populations).resolve().as_posix()} if populations else {}),
            "alignments" : [i.as_posix() for i in bamlist]
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column(header="value", justify="left")
    start_text.add_row("Samples:", f"{n}")
    start_text.add_row("Genome:", genome)
    start_text.add_row("Sample Pooling:", populations if populations else "no")
    start_text.add_row("Perform Phasing:", "yes" if vcf else "no")
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, "sv_naibr", start_text, output_dir, sm_log, quiet, "workflow/sv.naibr.summary")

sv.add_command(leviathan)
sv.add_command(naibr)
