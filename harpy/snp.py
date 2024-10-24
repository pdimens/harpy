"""Harpy workflows to call SNP variants"""

import os
import sys
import yaml
from pathlib import Path
from rich import box
from rich.table import Table
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, fetch_report, snakemake_log
from ._parsers import parse_alignment_inputs
from ._validations import check_fasta, validate_bam_RG, validate_popfile, validate_popsamples, validate_regions

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def snp():
    """
    Call SNPs and small indels on alignments
    
    Provide an additional subcommand `mpileup` or `freebayes` to get more information on using
    those variant callers. They are both robust variant callers, but `freebayes` is recommended when ploidy
    is greater than **2**.
    """

docstring = {
    "harpy snp mpileup": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--genome", "--ploidy", "--populations", "--regions"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ],
    "harpy snp freebayes": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--genome", "--ploidy", "--populations", "--regions"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/snp")
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False, readable=True), required = True, help = 'Genome assembly for variant calling')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "SNP/mpileup", show_default=True,  help = 'Output directory name')
@click.option('-n', '--ploidy', default = 2, show_default = True, type=click.IntRange(min = 1, max = 2), help = 'Ploidy of samples')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-r', '--regions', type=str, default=50000, show_default=True, help = "Regions where to call variants")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def mpileup(inputs, output_dir, regions, genome, threads, populations, ploidy, extra_params, snakemake, skip_reports, quiet, hpc, conda, setup_only):
    """
    Call variants from using bcftools mpileup
    
    Provide the input alignment (`.bam`) files and/or directories
    at the end of the command as individual files/folders, using shell wildcards
    (e.g. `data/scarab*.bam`), or both.
    
    The `--regions` option specifies what genomic regions to call variants
    with. If a BED or tab delimited file is provided, variant calling will be parallelized
    over those regions. If a single region is provided in the format `chrom:start-end`, only
    that region will be called. If an integer is provided (default), then Harpy will
    call variants in parallel for intervals of that size across the entire genome.

    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`. 
    """   
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/snp_mpileup.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_bam_RG(bamlist, threads, quiet)
    check_fasta(genome, quiet)

    # setup regions checks
    regtype = validate_regions(regions, genome)
    if regtype == "windows":
        region = Path(f"{workflowdir}/positions.bed").resolve()
        os.system(f"make_windows.py -m 1 -i {genome} -o {region} -w {regions}")
    elif regtype == "region":
        region = regions
    else:
        region = Path(f"{workflowdir}/positions.bed").resolve()
        os.system(f"cp -f {regions} {region}")

    fetch_rule(workflowdir, "snp_mpileup.smk")
    fetch_report(workflowdir, "bcftools_stats.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "snp_mpileup")
    if populations:
        validate_popfile(populations)
        # check that samplenames and populations line up
        validate_popsamples(bamlist, populations, quiet)
    configs = {
        "workflow" : "snp mpileup",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "ploidy" : ploidy,
        "region_type" : regtype,
        **({'windowsize': int(regions)} if regtype == "windows" else {}),
        **({'extra': extra_params} if extra_params else {}),
        "skip_reports" : skip_reports,
        "workflow_call" : command.rstrip(),
        "inputs" : {
            "genome" : Path(genome).resolve().as_posix(),
            "regions" : Path(region).resolve().as_posix() if regtype != "region" else region,
            **({'groupings': Path(populations).resolve().as_posix()} if populations else {}),
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
    start_text.add_column(header="value", justify="left")
    start_text.add_row("Samples:", f"{n}")
    if populations:
        start_text.add_row("Sample Groups:", populations)
    start_text.add_row("Genome:", genome)
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, "snp_mpileup", start_text, output_dir, sm_log, quiet, "workflow/snp.mpileup.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/snp")
@click.option('-x', '--extra-params', type = str, help = 'Additional variant caller parameters, in quotes')
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False, readable=True), required = True, help = 'Genome assembly for variant calling')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "SNP/freebayes", show_default=True,  help = 'Output directory name')
@click.option('-n', '--ploidy', default = 2, show_default = True, type=click.IntRange(min=1, max_open=True), help = 'Ploidy of samples')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True), help = "Tab-delimited file of sample\<tab\>population")
@click.option('-r', '--regions', type=str, default=50000, show_default=True, help = "Regions where to call variants")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def freebayes(inputs, output_dir, genome, threads, populations, ploidy, regions, extra_params, snakemake, skip_reports, quiet, hpc, conda, setup_only):
    """
    Call variants using freebayes
    
    Provide the input alignment (`.bam`) files and/or directories
    at the end of the command as individual files/folders, using shell wildcards
    (e.g. `data/jellyfish*.bam`), or both.
    
    The `--regions` option specifies what genomic regions to call variants
    with. If a BED or tab delimited file is provided, variant calling will be parallelized
    over those regions. If a single region is provided in the format `chrom:start-end`, only
    that region will be called. If an integer is provided (default), then Harpy will
    call variants in parallel for intervals of that size across the entire genome.

    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/snp_freebayes.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_bam_RG(bamlist, threads, quiet)
    check_fasta(genome, quiet)

    # setup regions checks
    regtype = validate_regions(regions, genome)
    if regtype == "windows":
        region = Path(f"{workflowdir}/positions.bed").resolve().as_posix()
        os.system(f"make_windows.py -m 0 -i {genome} -o {region} -w {regions}")
    elif regtype == "region":
        region = regions
    else:
        region = Path(f"{workflowdir}/positions.bed").resolve().as_posix()
        os.system(f"cp -f {regions} {region}")

    fetch_rule(workflowdir, "snp_freebayes.smk")
    fetch_report(workflowdir, "bcftools_stats.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "snp_freebayes")
    if populations:
        # check for delimeter and formatting
        validate_popfile(populations)
        # check that samplenames and populations line up
        validate_popsamples(bamlist, populations,quiet)

    configs = {
        "workflow" : "snp freebayes",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "ploidy" : ploidy,
        "region_type" : regtype,
        **({'windowsize': int(regions)} if regtype == "windows" else {}),
        **({'extra': extra_params} if extra_params else {}),
        "skip_reports" : skip_reports,
        "workflow_call" : command.rstrip(),
        "inputs" : {
            "genome" : Path(genome).resolve().as_posix(),
            "regions" : Path(region).resolve().as_posix() if regtype != "region" else region,
            **({'groupings': Path(populations).resolve().as_posix()} if populations else {}),
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
    start_text.add_column(header="value", justify="left")
    start_text.add_row("Samples:", f"{n}")
    if populations:
        start_text.add_row("Sample Groups:", populations)
    start_text.add_row("Genome:", genome)
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, "snp_freebayes", start_text, output_dir, sm_log, quiet, "workflow/snp.freebayes.summary")

snp.add_command(mpileup)
snp.add_command(freebayes)
