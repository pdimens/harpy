"""Harpy haplotype phasing workflow"""

import os
import sys
import yaml
from rich import box
from rich.table import Table
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, fetch_report, snakemake_log
from ._parsers import parse_alignment_inputs
from ._validations import check_fasta, vcf_samplematch, validate_bam_RG, validate_input_by_ext

docstring = {
        "harpy phase": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--genome", "--ignore-bx", "--molecule-distance", "--prune-threshold", "--vcf", "--vcf-samples"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },     
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/phase")
@click.option('-x', '--extra-params', type = str, help = 'Additional HapCut2 parameters, in quotes')
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False, readable=True), help = 'Path to genome assembly if wanting to also extract reads spanning indels')
@click.option('-b', '--ignore-bx',  is_flag = True, show_default = True, default = False, help = 'Ignore barcodes when phasing')
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Distance cutoff to split molecules (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Phase", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prune-threshold', default = 7, show_default = True, type = click.IntRange(0,100), help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 2, max_open = True), help = 'Number of threads to use')
@click.option('-v', '--vcf', required = True, type=click.Path(exists=True, dir_okay=False, readable=True), help = 'Path to BCF/VCF file')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake',  type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for phasing rather than those found the inputs')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def phase(inputs, output_dir, vcf, threads, molecule_distance, prune_threshold, vcf_samples, genome, snakemake, extra_params, ignore_bx, skip_reports, quiet, hpc, conda, setup_only):
    """
    Phase SNPs into haplotypes

    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/myotis*.bam`), or both.
    
    You may choose to omit barcode information with `--ignore-bx`, although it's usually
    better to include that information. Use `--vcf-samples` to phase only
    the samples present in your input `--vcf` file rather than all the samples present in
    the `INPUT` alignments.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/phase.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake:
        command += snakemake

    os.makedirs(f"{workflowdir}/input", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    samplenames = vcf_samplematch(vcf, bamlist, vcf_samples)
    validate_input_by_ext(vcf, "--vcf", ["vcf", "bcf", "vcf.gz"])
    validate_bam_RG(bamlist, threads, quiet)
    if genome:
        check_fasta(genome, quiet)
    fetch_rule(workflowdir, "phase.smk")
    fetch_report(workflowdir, "hapcut.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "phase")
    configs = {
        "workflow" : "phase",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "ignore_bx" : ignore_bx,
        "prune" : prune_threshold/100,
        "molecule_distance" : molecule_distance,
        "samples_from_vcf" : vcf_samples,
        **({'extra': extra_params} if extra_params else {}),
        "skip_reports" : skip_reports,
        "workflow_call" : command.rstrip(),
        "inputs" : {
            "variantfile" : vcf,
            **({'genome': genome} if genome else {}),
            "alignments" : [i.as_posix() for i in bamlist]
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes()
    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column("value", justify="left")
    start_text.add_row("Input VCF:", vcf)
    start_text.add_row("Samples in VCF:", f"{len(samplenames)}")
    start_text.add_row("Alignment Files:", f"{n}")
    start_text.add_row("Phase Indels:", "yes" if genome else "no")
    if genome:
        start_text.add_row("Genome:", genome)
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, "phase", start_text, output_dir, sm_log, quiet, "workflow/phase.summary")
