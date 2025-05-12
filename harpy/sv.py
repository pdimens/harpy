"""Harpy workflows to detect structural variants"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._cli_types_generic import convert_to_int, ContigList, HPCProfile, InputFile, SnakemakeParams
from ._cli_types_params import LeviathanParams, NaibrParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, fetch_report, instantiate_dir, setup_snakemake, write_workflow_config
from ._parsers import parse_alignment_inputs
from ._printing import workflow_info
from ._validations import check_fasta, check_phase_vcf
from ._validations import validate_popfile, validate_popsamples, fasta_contig_match

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def sv():
    """
    Call large structural variants on alignments
 
    | caller | inversions | duplications | deletions | breakends |
    |:-------|:----------:|:------------:|:---------:|:---------:|
    | leviathan |      ✔  |     ✔        |     ✔     |      ✔    |
    | naibr     |      ✔  |     ✔        |     ✔     |     🗙    |

    Provide the subcommand `leviathan` or `naibr` to get more information on using
    those variant callers. NAIBR tends to call variants better, but requires more user preprocessing.
    """

module_docstring = {
    "harpy sv": [
        {
            "name": "Commands",
            "commands": ["leviathan", "naibr"],
            "panel_styles": {"border_style": "blue"}
        }
    ]
}

docstring = {
    "harpy sv leviathan": [
        {
            "name": "Parameters",
            "options": ["--duplicates", "--extra-params", "--iterations", "--min-barcodes", "--min-size", "--populations", "--sharing-thresholds"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy sv naibr": [
        {
            "name": "Module Parameters",
            "options": ["--extra-params", "--min-barcodes", "--min-quality", "--min-size", "--molecule-distance", "--populations", "--vcf"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
} | module_docstring

@click.command(no_args_is_help = True, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/sv/leviathan/")
@click.option('-x', '--extra-params', type = LeviathanParams(), help = 'Additional leviathan parameters, in quotes')
@click.option('-i', '--iterations', show_default = True, default=50, type = click.IntRange(min = 10), help = 'Number of iterations to perform through index (reduces memory)')
@click.option('-d', '--duplicates', show_default = True, default=10, type = click.IntRange(min = 1), help = 'Consider SV of the same type as duplicates if their breakpoints are within this distance')
@click.option('-m', '--min-size', type = click.IntRange(min = 10), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-s', '--sharing-thresholds', type = click.IntRange(0,100, clamp = True), default = (95,95,95), nargs = 3, show_default=True, help = 'Percentile thresholds in the distributions of the number of shared barcodes for (small,medium,large) variants, separated by spaces')
@click.option('-b', '--min-barcodes', show_default = True, default=2, type = click.IntRange(min = 1), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "SV/leviathan", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`\\<TAB\\>`population`')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=InputFile("fasta", gzip_ok = True), required = True, nargs = 1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def leviathan(inputs, output_dir, reference, min_size, min_barcodes, iterations, duplicates, sharing_thresholds, threads, populations, extra_params, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Call structural variants using LEVIATHAN
    
    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    Optionally specify `--populations` for population-pooled variant calling
    (**harpy template** can create that file). If you suspect Leviathan is missing certain variants
    you expect to find, try lowering `--sharing-thresholds`, _e.g._ `95,95,95`. The thresholds don't
    have to be the same across the different size classes.
    """
    workflow = "sv_leviathan"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    bamlist, n = parse_alignment_inputs(inputs, "INPUTS")
    check_fasta(reference)
    if contigs:
        fasta_contig_match(contigs, reference)
    if populations:
        validate_popfile(populations)
        validate_popsamples(bamlist, populations,quiet)

    ## setup workflow ##
    vcaller = workflow if not populations else f"{workflow}_pop"
    command,command_rel = setup_snakemake(
        vcaller,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, f"{vcaller}.smk")
    fetch_report(workflowdir, "leviathan.qmd")
    fetch_report(workflowdir, "leviathan_pop.qmd") if populations else None

    conda_envs = ["align", "r", "variants"]
    configs = {
        "workflow" : workflow,
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
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "reference" : reference,
            **({'groupings': populations} if populations else {}),
            "alignments" : bamlist
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:", n),
        ("Reference:", os.path.basename(reference)),
        ("Sample Pooling:", os.path.basename(populations) if populations else "no"),
        ("Output Folder:", os.path.basename(output_dir) + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/sv.leviathan.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/sv/naibr/")
@click.option('-x', '--extra-params', type = NaibrParams(), help = 'Additional naibr parameters, in quotes')
@click.option('-b', '--min-barcodes', show_default = True, default=2, type = click.IntRange(min = 1), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-q', '--min-quality', show_default = True, default=30, type = click.IntRange(min = 0, max = 40), help = 'Minimum mapping quality of reads to use')
@click.option('-m', '--min-size', type = click.IntRange(min = 10), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Base-pair distance delineating separate molecules')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "SV/naibr", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`\\<TAB\\>`population`')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-v', '--vcf', type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True),  help = 'Path to phased bcf/vcf file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=InputFile("fasta", gzip_ok = True), required = True, nargs = 1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def naibr(inputs, output_dir, reference, vcf, min_size, min_barcodes, min_quality, threads, populations, molecule_distance, extra_params, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Call structural variants using NAIBR
    
    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    NAIBR requires **phased** bam files as input. This appears as the `HP` or `PS` tags
    in alignment records. If your bam files do not have either of these phasing tags
    (e.g. BWA/EMA do not phase alignments), then provide a **phased** `--vcf` file such
     as that created by `harpy phase` and Harpy will use [whatshap haplotag](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)
    to phase your input bam files prior to calling variants with NAIBR.

    Optionally specify `--populations` for population-pooled variant calling (**harpy template** can create that file).
    """
    workflow = "sv_naibr"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    bamlist, n = parse_alignment_inputs(inputs, "INPUTS")
    check_fasta(reference)
    if contigs:
        fasta_contig_match(contigs, reference)
    if populations:
        validate_popfile(populations)
        validate_popsamples(bamlist, populations, quiet)
    if vcf:
        check_phase_vcf(vcf)

    ## setup workflow ##
    vcaller = workflow if not populations else f"{workflow}_pop"
    vcaller += "_phase" if vcf else ""
    command,command_rel = setup_snakemake(
        vcaller,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, f"{vcaller}.smk")
    fetch_report(workflowdir, "naibr_pop.qmd") if populations else None
    fetch_report(workflowdir, "naibr.qmd")

    conda_envs = ["phase", "r", "variants"]
    configs = {
        "workflow" : workflow,
        "min_barcodes" : min_barcodes,
        "min_quality" : min_quality,
        "min_size" : min_size,
        "molecule_distance" : molecule_distance,
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            **({'reference': reference} if reference else {}),
            **({'vcf': vcf} if vcf else {}),
            **({'groupings': populations} if populations else {}),
            "alignments" : bamlist
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Samples:", n),
        ("Reference:", os.path.basename(reference)),
        ("Sample Pooling:", os.path.basename(populations) if populations else "no"),
        ("Perform Phasing:", "yes" if vcf else "no"),
        ("Output Folder:", os.path.basename(output_dir) + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/sv.naibr.summary")

sv.add_command(leviathan)
sv.add_command(naibr)
