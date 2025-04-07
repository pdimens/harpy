"""Harpy workflows to simulate genomic variants and linked-reads"""
import os
import sys
import yaml
import shutil
from pathlib import Path
import rich_click as click
from ._cli_types_generic import convert_to_int, HPCProfile, InputFile, SnakemakeParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, fetch_script, snakemake_log, write_snakemake_config, write_workflow_config
from ._printing import workflow_info
from ._validations import check_fasta, validate_barcodefile

docstring = {
    "harpy simulate linkedreads": [
        {
            "name": "Parameters",
            "options": ["--barcodes", "--distance-sd", "--outer-distance", "--molecule-length", "--molecules-per", "--mutation-rate", "--partitions", "--read-pairs"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--merge-haplotypes", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        }
    ]
}

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/simulate/simulate-linkedreads")
@click.option('-b', '--barcodes', type = click.Path(exists=True, dir_okay=False, readable=True), help = "File of linked-read barcodes to add to reads")
@click.option('-s', '--distance-sd', type = click.IntRange(min = 1), default = 15, show_default=True,  help = "Standard deviation of read-pair distance")
@click.option('-m', '--molecules-per', type = click.IntRange(min = 1, max = 4700), default = 10, show_default=True,  help = "Average number of molecules per partition")
@click.option('-l', '--molecule-length', type = click.IntRange(min = 2), default = 100, show_default=True,  help = "Mean molecule length (kbp)")
@click.option('-r', '--mutation-rate', type = click.FloatRange(min = 0), default=0.001, show_default=True,  help = "Random mutation rate for simulating reads")
@click.option('-d', '--outer-distance', type = click.IntRange(min = 100), default = 350, show_default= True, help = "Outer distance between paired-end reads (bp)")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Simulate/linkedreads", help = 'Output directory name')
@click.option('-p', '--partitions', type = click.IntRange(min = 1), default=1500, show_default=True,  help = "Number (in thousands) of partitions/beads to generate")
@click.option('-n', '--read-pairs', type = click.FloatRange(min = 0.001), default = 600, show_default=True,  help = "Number (in millions) of read pairs to simulate")
@click.option('--merge-haplotypes',  is_flag = True, default = False, help = 'Concatenate sequences across haplotypes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('genome_hap1', required=True, type = InputFile("fasta", gzip_ok = True), nargs=1)
@click.argument('genome_hap2', required=True, type = InputFile("fasta", gzip_ok = True), nargs=1)
def linkedreads(genome_hap1, genome_hap2, output_dir, outer_distance, mutation_rate, distance_sd, barcodes, read_pairs, molecule_length, partitions, molecules_per, threads, snakemake, quiet, merge_haplotypes, hpc, container, setup_only):
    """
    Create linked reads from a genome
 
    Two haplotype genomes (un/compressed fasta) need to be provided as inputs at the end of the command. If
    you don't have a diploid genome, you can simulate one with `harpy simulate` as described [in the documentation](https://pdimens.github.io/harpy/blog/simulate_diploid/).

    If not providing a file for `--barcodes`, Harpy will generate a file containing the original
    (96^4) set of 24-basepair haplotagging barcodes (~2GB disk space). The `--barcodes` file is expected to have one
    linked-read barcode per line, given as nucleotides. Use `--merge-haplotypes` to merge haplotype 1 and haplotype 2 for R1 reads
    (same for R2), resulting in one file of R1 reads and one file of R2 reads.
    """
    ## checks and validations ##
    check_fasta(genome_hap1)
    check_fasta(genome_hap2)
    if barcodes:
        bc_len = validate_barcodefile(barcodes, True, quiet, gzip_ok=False, haplotag_only=True)

    ## setup workflow ##
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    write_snakemake_config("conda" if not container else "conda apptainer", output_dir)
    command = f"snakemake --cores {threads} --snakefile {workflowdir}/simulate_linkedreads.smk"
    command += f" --configfile {workflowdir}/config.harpy.yaml --profile {workflowdir}"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    fetch_rule(workflowdir, "simulate_linkedreads.smk")
    fetch_script(workflowdir, "HaploSim.pl")

    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "simulate_linkedreads")

    conda_envs = ["simulations"]
    configs = {
        "workflow" : "simulate linkedreads",
        "snakemake_log" : sm_log,
        "outer_distance" : outer_distance,
        "distance_sd" : distance_sd,
        "read_pairs" : read_pairs,
        "mutation_rate" : mutation_rate,
        "molecule_length" : molecule_length,
        "partitions" : partitions,
        "molecules_per_partition" : molecules_per,
        "merge_haplotypes": merge_haplotypes,
        "snakemake_command" : command.rstrip(),
        "conda_environments" : conda_envs,
        'barcodes': {
            "file": Path(barcodes).resolve().as_posix() if barcodes else f"{workflowdir}/input/haplotag_barcodes.txt",
            "length": bc_len if barcodes else 24
        },
        "inputs" : {
            "genome_hap1" : Path(genome_hap1).resolve().as_posix(),
            "genome_hap2" : Path(genome_hap2).resolve().as_posix(),
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Genome Haplotype 1:", os.path.basename(genome_hap1)),
        ("Genome Haplotype 2:", os.path.basename(genome_hap2)),
        ("Barcodes:", os.path.basename(barcodes) if barcodes else "Haplotagging Default"),
        ("Output Folder:", output_dir + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "simulate_linkedreads", start_text, output_dir, sm_log, quiet, "workflow/simulate.reads.summary")

