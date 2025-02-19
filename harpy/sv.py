"""Harpy workflows to detect structural variants"""

import os
import sys
import yaml
from pathlib import Path
import rich_click as click
from ._cli_types_generic import convert_to_int, ContigList, HPCProfile, InputFile, IntList, SnakemakeParams
from ._cli_types_params import LeviathanParams, NaibrParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake, SNAKEMAKE_CMD
from ._misc import fetch_rule, fetch_report, snakemake_log
from ._parsers import parse_alignment_inputs
from ._printing import workflow_info
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
            "options": ["--duplicates", "--extra-params", "--genome", "--iterations", "--min-barcodes", "--min-size", "--populations", "--sharing-thresholds"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--container", "--contigs", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ],
    "harpy sv naibr": [
        {
            "name": "Module Parameters",
            "options": ["--extra-params", "--genome", "--min-barcodes", "--min-quality", "--min-size", "--molecule-distance", "--populations", "--vcf"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--container", "--contigs", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/sv/leviathan/")
@click.option('-x', '--extra-params', type = str, help = 'Additional leviathan parameters, in quotes')
@click.option('-g', '--genome', type=InputFile("fasta", gzip_ok = True), required = True, help = 'Genome assembly for variant calling')
@click.option('-i', '--iterations', show_default = True, default=50, type = click.IntRange(min = 10, max_open = True), help = 'Number of iterations to perform through index (reduces memory)')
@click.option('-d', '--duplicates', show_default = True, default=10, type = click.IntRange(min = 1, max_open = True), help = 'Consider SV of the same type as duplicates if their breakpoints are within this distance')
@click.option('-s', '--min-size', type = click.IntRange(min = 10, max_open = True), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-r', '--sharing-thresholds', type = IntList(3), default = "99,99,99", show_default=True, help = 'Percentile thresholds in the distributions of the number of shared barcodes for (small,medium,large) variants, separated by commas')
@click.option('-b', '--min-barcodes', show_default = True, default=2, type = click.IntRange(min = 1, max_open = True), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "SV/leviathan", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True), help = 'File of `sample`\\<TAB\\>`population`')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max = 999), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def leviathan(inputs, output_dir, genome, min_size, min_barcodes, iterations, duplicates, sharing_thresholds, threads, populations, extra_params, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Call structural variants using LEVIATHAN
    
    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    Optionally specify `--populations` for population-pooled variant calling
    (**harpy popgroup** can create that file). If you suspect Leviathan is missing certain variants
    you expect to find, try lowering `--sharing-thresholds`, _e.g._ `95,95,95`. The thresholds don't
    have to be the same across the different size classes.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    vcaller = "leviathan" if populations is None else "leviathan_pop"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/sv_{vcaller}.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        command += f" --workflow-profile {hpc}"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    check_fasta(genome)
    if contigs:
        fasta_contig_match(contigs, genome)
    if populations:
        validate_popfile(populations)
        validate_popsamples(bamlist, populations,quiet)
        fetch_report(workflowdir, "leviathan_pop.qmd")
    fetch_report(workflowdir, "leviathan.qmd")
    fetch_rule(workflowdir, f"sv_{vcaller}.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "sv_leviathan")
    conda_envs = ["align", "r", "variants"]
    configs = {
        "workflow" : "sv leviathan",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "min_barcodes" : min_barcodes,
        "min_size" : min_size,
        "iterations" : iterations,
        "variant_thresholds": {
            "small" : sharing_thresholds[0],
            "medium" : sharing_thresholds[1],
            "large" : sharing_thresholds[2],
            "duplicates": duplicates
        },
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

    start_text = workflow_info(
        ("Samples:", n),
        ("Genome:", genome),
        ("Sample Pooling:", populations if populations else "no"),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "sv_leviathan", start_text, output_dir, sm_log, quiet, "workflow/sv.leviathan.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/sv/naibr/")
@click.option('-x', '--extra-params', type = str, help = 'Additional naibr parameters, in quotes')
@click.option('-g', '--genome', type=InputFile("fasta", gzip_ok = True), required = True, help = 'Genome assembly for calling variants')
@click.option('-b', '--min-barcodes', show_default = True, default=2, type = click.IntRange(min = 1, max_open = True), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-q', '--min-quality', show_default = True, default=30, type = click.IntRange(min = 0, max = 40), help = 'Minimum mapping quality of reads to use')
@click.option('-s', '--min-size', type = click.IntRange(min = 10, max_open = True), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Base-pair distance delineating separate molecules')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "SV/naibr", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True), help = 'File of `sample`\\<TAB\\>`population`')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max = 999), help = 'Number of threads to use')
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True),  help = 'Path to phased bcf/vcf file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def naibr(inputs, output_dir, genome, vcf, min_size, min_barcodes, min_quality, threads, populations, molecule_distance, extra_params, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Call structural variants using NAIBR
    
    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    NAIBR requires **phased** bam files as input. This appears as the `HP` or `PS` tags
    in alignment records. If your bam files do not have either of these phasing tags
    (e.g. BWA/EMA do not phase alignments), then provide a **phased** `--vcf` file such
     as that created by `harpy phase` and Harpy will use [whatshap haplotag](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)
    to phase your input bam files prior to calling variants with NAIBR.

    Optionally specify `--populations` for population-pooled variant calling (**harpy popgroup** can create that file).
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    vcaller = "naibr" if populations is None else "naibr_pop"
    vcaller += "_phase" if vcf else ""
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/sv_{vcaller}.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        command += f" --workflow-profile {hpc}"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    check_fasta(genome)
    if contigs:
        fasta_contig_match(contigs, genome)
    if populations:
        validate_popfile(populations)
        validate_popsamples(bamlist, populations, quiet)
        fetch_report(workflowdir, "naibr_pop.qmd")
    fetch_report(workflowdir, "naibr.qmd")
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
        "min_size" : min_size,
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

    start_text = workflow_info(
        ("Samples:", n),
        ("Genome:", genome),
        ("Sample Pooling:", populations if populations else "no"),
        ("Perform Phasing:", "yes" if vcf else "no"),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "sv_naibr", start_text, output_dir, sm_log, quiet, "workflow/sv.naibr.summary")

sv.add_command(leviathan)
sv.add_command(naibr)
