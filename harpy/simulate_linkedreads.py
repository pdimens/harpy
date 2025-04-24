"""Harpy workflows to simulate genomic variants and linked-reads"""
import os
import sys
import yaml
import shutil
import rich_click as click
from ._cli_types_generic import convert_to_int, HPCProfile, InputFile, SnakemakeParams
from ._cli_types_params import Barcodes
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, instantiate_dir, setup_snakemake, write_workflow_config
from ._printing import workflow_info
from ._validations import check_fasta

docstring = {
    "harpy simulate linkedreads": [
        {
            "name": "Read Simulation Parameters",
            "options": ["--coverage","--distance","--error", "--length",  "--regions", "--stdev"],
            "panel_styles": {"border_style": "dim blue"}
        },
        {
            "name": "Linked Read Parameters",
            "options": ["--lr-type", "--molecule-coverage", "--molecule-length", "--molecule-number"],
            "panel_styles": {"border_style": "dim magenta"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--help", "--hpc", "--output-prefix", "--output-type", "--quiet", "--snakemake", "--threads", "--version"],
            "panel_styles": {"border_style": "dim"}
        }
    ]
}

@click.command(epilog = "Documentation: https://pdimens.github.io/mimick/", no_args_is_help = True)
#Paired-end FASTQ simulation using pywgsim
@click.option('--coverage', help='mean coverage target for simulated data', show_default=True, default=30, type=click.FloatRange(min=0.05))
@click.option('--distance', help='outer distance between the two ends in bp', default=500, show_default=True, type=click.IntRange(min=0))
@click.option('--error', help='base error rate', default=0.02, show_default=True, type=click.FloatRange(min=0, max=1))
@click.option('--extindels', help='indels extension rate', default=0, hidden=True, type=click.FloatRange(min=0, max=1))
@click.option('--indels', help='indels rate', default=0, hidden=True, type=click.FloatRange(min=0, max=1))
@click.option('--length', help='length of reads in bp', default=150, show_default=True, type=click.IntRange(min=30))
@click.option('--mutation', help='mutation rate', default=0, hidden=True, type=click.FloatRange(min=0))
@click.option('--regions', help='one or more regions to simulate, in BED format', type = click.Path(dir_okay=False, readable=True, resolve_path=True))
@click.option('--stdev', help='standard deviation of --distance', default=50, show_default=True, type=click.IntRange(min=0))
#Linked-read simulation
@click.option('-l', '--lr-type', help='type of linked-read experiment', default = "haplotagging", show_default=True, show_choices=True, type= click.Choice(["10x", "stlfr", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-c','--molecule-coverage', help='mean percent coverage per molecule if <1, else mean number of reads per molecule', default=0.2, show_default=True, type=click.FloatRange(min=0.00001))
@click.option('-m','--molecule-length', help='mean length of molecules in bp', show_default=True, default=80000, type=click.IntRange(min=50))
@click.option('-n','--molecule-number', help='mean number of unrelated molecules per barcode, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a Poisson distribution', default=3, show_default=True, type=int)
# general
@click.option('-o','--output-prefix', help='output file prefix', type = click.Path(exists = False, writable=True, resolve_path=True), default = "Simulate/linkedreads/SIM", show_default=True)
@click.option('-O','--output-type', help='output format of FASTQ files', type = click.Choice(["10x", "stlfr", "standard", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-t','--threads', help='number of threads to use for simulation', type = click.IntRange(1, 999, clamp = True), default=2, show_default=True)
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('barcodes', type = Barcodes())
@click.argument('fasta', type = InputFile("fasta", gzip_ok = True), nargs = -1, required=True)
def linkedreads(barcodes, fasta, output_prefix, output_type, regions, threads,coverage,distance,error,extindels,indels,length,mutation,stdev,lr_type, molecule_coverage, molecule_length, molecule_number, snakemake, quiet, hpc, container, setup_only):
    """
    Create linked reads using genome haplotypes

    Barcodes can be supplied one of two ways:
   
    1. randomly generate barcodes based on a specification of `length,count`
        - two integers, comma-separated, no space
        - e.g. `16,400000` would generate 400,000 unique 16bp barcodes 
    2. providing a file of nucleotide barcodes, 1 per line

    You can specify the linked-read barcode chemistry to simulate via `--lr-type` as well as
    the output format of FASTQ files (default: same as `lr-type`).

    | --lr-type | Format |
    |:------------------|:-------|
    |`10x`/`tellseq`   | single barcode on R1 |
    |`haplotagging`  | R1 and R2 each have different combinatorial 2-barcodes |
    |`stlfr`         | combinatorial 3-barcode on R2 |

    | --output-type | Barcode Location | Example |
    |:-----------------|:-------|:---------------------|
    |`10x`           | start of R1 sequence | `ATAGACCATAGA`GGACA... |
    |`haplotagging`  | sequence header as `BX:Z:ACBD` |  `@SEQID BX:Z:A0C331B34D87` |
    |`standard`      | sequence header as `BX:Z:BARCODE`, no specific format | `@SEQID BX:Z:ATACGAGACA` |
    |`stlfr`         | appended to sequence ID via `#1_2_3` | `@SEQID#1_354_39` |
    |`tellseq`       | appended to sequence ID via `:ATCG` | `@SEQID:TATTAGCAC` |
    """
    workflow = "simulate_linkedreads"
    output_dir = os.path.dirname(output_prefix)
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    for i in fasta:
        check_fasta(i)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        workflow,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "simulate_linkedreads.smk")

    conda_envs = ["simulations"]
    configs = {
        "workflow" : workflow,
        "read_coverage" : coverage,
        "outer_distance" : distance,
        "error_rate" :   error,
        "length" :  length,
        "stdev" : stdev,
        "lr-type" : lr_type,            
        "molecule-coverage" : molecule_coverage,  
        "molecule-length" : molecule_length,
        "molecule-number" : molecule_number,   
        "mutation" : mutation,
        "indels" :  indels,
        "extindels" : extindels,
        "output-prefix" : output_prefix,
        "output-type" : output_type if output_type else lr_type,
        **({"regions":  regions} if regions else {}),
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "inputs" : {
            "barcodes" : barcodes,
            "haplotypes" : list(fasta),
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Haplotypes:", len(fasta)),
        ("Barcodes:", os.path.basename(barcodes) if os.path.exists(barcodes) else "randomly generated"),
        ("Output Folder:", os.path.basename(output_dir) + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/simulate.reads.summary")
