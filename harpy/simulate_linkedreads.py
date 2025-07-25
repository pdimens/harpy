"""Harpy workflows to simulate genomic variants and linked reads"""
import os
import sys
import yaml
import shutil
import rich_click as click
from ._cli_types_generic import HPCProfile, InputFile, ReadLengths, SnakemakeParams
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
            "options": ["--coverage","--distance","--error", "--lengths",  "--regions", "--stdev"],
            "panel_styles": {"border_style": "dim blue"}
        },
        {
            "name": "Linked Read Parameters",
            "options": ["--circular", "--segments", "--molecule-attempts", "--molecule-coverage", "--molecule-length", "--molecules-per", "--singletons"],
            "panel_styles": {"border_style": "dim magenta"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--help", "--hpc", "--output-prefix", "--output-type", "--quiet", "--seed", "--snakemake", "--threads", "--version"],
            "panel_styles": {"border_style": "dim"}
        }
    ]
}

@click.command(epilog = "Documentation: https://pdimens.github.io/mimick/", no_args_is_help = True)
#Paired-end FASTQ simulation using pywgsim
@click.option('--coverage', help='mean coverage (depth) target for simulated data', show_default=True, default=30, type=click.FloatRange(min=0.05))
@click.option('--distance', help='outer distance between the two read ends in bp', default=500, show_default=True, type=click.IntRange(min=0))
@click.option('--error', help='base error rate', default=0.02, show_default=True, type=click.FloatRange(min=0, max=1))
@click.option('--extindels', help='indels extension rate', default=0, hidden=True, type=click.FloatRange(min=0, max=1))
@click.option('--indels', help='indels rate', default=0, hidden=True, type=click.FloatRange(min=0, max=1))
@click.option('--lengths', help='length of R1,R2 reads in bp', default="150,150", show_default=True, type=ReadLengths())
@click.option('--mutation', help='mutation rate', default=0, hidden=True, type=click.FloatRange(min=0))
@click.option('--regions', help='one or more regions to simulate, in BED format', type = click.Path(dir_okay=False, readable=True, resolve_path=True))
@click.option('--stdev', help='standard deviation for `--distance`', default=50, show_default=True, type=click.IntRange(min=0))
#Linked-read simulation
@click.option('-C','--circular', is_flag= True, default = False, help = 'contigs are circular/prokaryotic')
@click.option('-x', '--segments', help='treat barcodes as combinatorial with this many segments', default = 4, show_default=True, type= click.IntRange(min=1))
@click.option('-a','--molecule-attempts', help='how many tries to create a molecule with <70% ambiguous bases', show_default=True, default=300, type=click.IntRange(min=5))
@click.option('-c','--molecule-coverage', help='mean percent coverage per molecule if <1, else mean number of reads per molecule', default=0.2, show_default=True, type=click.FloatRange(min=0.00001))
@click.option('-m','--molecule-length', help='mean length of molecules in bp', show_default=True, default=80000, type=click.IntRange(min=50))
@click.option('-n','--molecules-per', help='mean number of unrelated molecules per barcode, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a Normal distribution', default=2, show_default=True, type=int)
@click.option('-s','--singletons', help='proportion of barcodes that will only have one read pair', default=0, show_default=True, type=click.FloatRange(0,1))
# general
@click.option('-o','--output-prefix', help='output file prefix', type = click.Path(exists = False, writable=True, resolve_path=True), default = "Simulate/linkedreads/SIM", show_default=True)
@click.option('-O','--output-type', help='FASTQ output format', type = click.Choice(["10x", "stlfr", "standard", "standard:haplotagging", "standard:stlfr", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-S','--seed', help='random seed for simulation', type=click.IntRange(min=1), default=None)
@click.option('-t','--threads', help='number of threads to use for simulation', type = click.IntRange(1, 999, clamp = True), default=2, show_default=True)
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('barcodes', type = Barcodes())
@click.argument('fasta', type = InputFile("fasta", gzip_ok = True), nargs = -1, required=True)
def linkedreads(barcodes, fasta, output_prefix, output_type, regions, threads,coverage,distance,error,extindels,indels,lengths,mutation,stdev,circular,segments, molecule_attempts,molecule_coverage, molecule_length, molecules_per, singletons, seed, snakemake, quiet, hpc, container, setup_only):
    """
    Create linked reads using genome haplotypes

    This workflow is a very thin veneer to use `Mimick` with mutations (SNPs and indels) deactivated.
    For all the features, you are encouraged to install `Mimick` from Bioconda and use it directly.
    In addition to selecting an `--output-type` (default varies by `-x`), barcodes can be parsed absolutely or you can specify the
    linked-read barcode type using `-x/--segments`. For example, to simulate the common 4-segment haplotagging style:
    ```
    harpy simulate linkedreads -x 4 -O haplotagging 6,96 hap.{1..2}.fa
    ```
    
    The `standard` output type can be suffixed with `:haplotagging` or `:stlfr` to use those barcode styles with the standard format
    (e.g. `standard:haplotagging`). See the [Mimick documentation](https://pdimens.github.io/mimick/#/data_formats)
    for a thorough explanation of the chemistries and output formats.
    Configurations for the common linked-read varieties: 

    | chemistry    | `--segments` | `--lengths` | `--output-type` default |
    |:-------------|:------------:|:-----------:|:------------------------|
    | 10x          |     `1`      |  `134,150`  | `tellseq`               |
    | tellseq      |     `1`      |  `132,150`  | `tellseq`               |
    | haplotagging |     `4`      |  `150,150`  | `standard:haplotagging` |
    | stlfr        |     `3`      |  `150,108`  | `stlfr`                 |
   
    Input barcodes can be supplied one of two ways:
    1. randomly generate barcodes based on a specification of `length,count`
        - two integers, comma-separated, no space
        - e.g. `16,400000` would generate 400,000 unique 16bp barcodes 
    2. providing a file of nucleotide barcodes, 1 per line
        ```
        ATGCAGGA
        GGAGGACT
        ```

    """
    workflow = "simulate_linkedreads"
    output_dir = os.path.dirname(output_prefix)
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    for i in fasta:
        check_fasta(i)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
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
        "output-prefix" : os.path.basename(output_prefix),
        "read_params": {
            "read_coverage" : coverage,
            "outer_distance" : distance,
            "error_rate" :   error,
            "lengths" :  lengths,
            "stdev" : stdev
        },
        "variant_params":{
            "mutation" : mutation,
            "indels" :  indels,
            "extindels" : extindels
        },
        "linked_read_params": {
            "circular": circular,
            "segments" : segments,            
            "output-type" : output_type if output_type else None,
            "molecule-attempts": molecule_attempts,
            "molecule-coverage" : molecule_coverage,  
            "molecule-length" : molecule_length,
            "molecules-per" : molecules_per,
            "singletons": singletons,
        },
        "random_seed": seed,
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
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/simulate.reads.summary")
