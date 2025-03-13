"""Harpy haplotype phasing workflow"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake, SNAKEMAKE_CMD
from ._misc import fetch_rule, fetch_report, snakemake_log
from ._cli_types_generic import convert_to_int, ContigList, HPCProfile, InputFile, SnakemakeParams
from ._cli_types_params import HapCutParams
from ._parsers import parse_alignment_inputs
from ._printing import workflow_info
from ._validations import check_fasta, vcf_sample_match, validate_bam_RG, vcf_contig_match

docstring = {
        "harpy phase": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--genome", "--ignore-bx", "--molecule-distance", "--prune-threshold", "--vcf", "--vcf-samples"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },     
    ]
}

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/phase")
@click.option('-x', '--extra-params', type = HapCutParams(), help = 'Additional HapCut2 parameters, in quotes')
@click.option('-g', '--genome', type=InputFile("fasta", gzip_ok = True), help = 'Path to genome assembly if wanting to also extract reads spanning indels')
@click.option('-b', '--ignore-bx',  is_flag = True, show_default = True, default = False, help = 'Ignore barcodes when phasing')
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = int, help = 'Distance cutoff to split molecules (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Phase", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prune-threshold', default = 7, show_default = True, type = click.IntRange(0,100, clamp = True), help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(2, 999, clamp = True), help = 'Number of threads to use')
@click.option('-v', '--vcf', required = True, type = InputFile("vcf", gzip_ok = False), help = 'Path to BCF/VCF file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for phasing rather than those found the inputs')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def phase(inputs, output_dir, vcf, threads, molecule_distance, prune_threshold, vcf_samples, genome, snakemake, extra_params, ignore_bx, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Configure and run haplotype phasing workflow using Snakemake.
    
    This function sets up the necessary workflow environment by creating required
    directories and generating a YAML configuration file based on provided inputs.
    It processes alignment files (BAMs) and a VCF file with variant calls while
    applying optional filters such as sample restriction and barcode omission.
    When an HPC configuration file is provided, it is copied into a dedicated subdirectory
    and referenced in the workflow command. If not in setup-only mode, the function
    launches the Snakemake workflow to phase SNPs into haplotypes.
    
    Args:
        inputs: Input BAM files or directories; supports shell wildcards.
        output_dir: Root directory for workflow outputs, configuration, and logs.
        vcf: VCF file containing variant calls used for phasing.
        threads: Number of CPU cores allocated for the workflow.
        molecule_distance: Maximum allowed distance between molecules for linking.
        prune_threshold: Pruning threshold as a percentage (divided by 100 in config).
        vcf_samples: If True, phase only the samples present in the VCF file.
        genome: Path to the genome FASTA file, necessary for indel phasing.
        snakemake: Additional command-line arguments to append to the Snakemake command.
        extra_params: Extra configuration parameters to merge into the YAML configuration.
        ignore_bx: If True, omit barcode information from the analysis.
        skip_reports: If True, skip generation of workflow reports.
        quiet: If True, reduce output verbosity.
        hpc: Path to an HPC submission YAML configuration file for remote execution.
        container: If True, use containerized execution (e.g., Apptainer) in the workflow.
        contigs: Optional contig specification for filtering or plotting VCF data.
        setup_only: If True, set up the workflow without launching the execution.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/phase.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    os.makedirs(f"{workflowdir}/input", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    samplenames = vcf_sample_match(vcf, bamlist, vcf_samples)
    validate_bam_RG(bamlist, threads, quiet)
    if genome:
        check_fasta(genome)
    if contigs:
        vcf_contig_match(contigs, vcf)
    fetch_rule(workflowdir, "phase.smk")
    fetch_report(workflowdir, "hapcut.qmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "phase")
    conda_envs = ["phase", "r"]
    configs = {
        "workflow" : "phase",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "ignore_bx" : ignore_bx,
        "prune" : prune_threshold/100,
        "molecule_distance" : molecule_distance,
        "samples_from_vcf" : vcf_samples,
        **({'extra': extra_params} if extra_params else {}),
        "workflow_call" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "variantfile" : vcf,
            **({'genome': genome} if genome else {}),
            "alignments" : [i.as_posix() for i in bamlist]
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Input VCF:", vcf),
        ("Samples in VCF:", len(samplenames)),
        ("Alignment Files:", n),
        ("Phase Indels:", "yes" if genome else "no"),
        ("Genome:", genome) if genome else None,
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "phase", start_text, output_dir, sm_log, quiet, "workflow/phase.summary")
