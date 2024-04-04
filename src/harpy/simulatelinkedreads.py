from .helperfunctions import generate_conda_deps, fetch_file 
from .printfunctions import print_onstart
import rich_click as click
import subprocess
import os
import sys
"""
    Reference genome and variants:
    -i INT      Outer distance between the two ends for pairs [350] XX
    -s INT      Standard deviation of the distance for pairs [35] XX
    -b STRING   Barcode list                                      XX
    -x INT      # million reads pairs in total to simulated [600] XX
    -f INT      Mean molecule length in kbp [100]                 XX
    -t INT      n*1000 partitions to generate [1500]              XX
    -m INT      Average # of molecules per partition [10]         XX
    
    -1 INT      1 SNP per INT base pairs [1000]
    -z INT      # of threads to run DWGSIM [8]
    -g STRING   Haploid FASTAs separated by comma. Overrides -r and -d.
    -d INT      Haplotypes to simulate [2]
    -o          Disable parameter checking
    -u 2
"""

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/simulate")
@click.option('-d', '--outer-distance', type = int, default = 350, show_default= True, help = "Outer distance between paired-end reads (bp)")
@click.option('-i', '--distance-sd', type = int, default = 15, show_default=True, help = "Standard deviation of read-pair distance")
@click.option('-b', '--barcodes', type = click.Path(exists=True), help = "File of linked-read barcodes")
@click.option('-n', '--read-pairs', type = int, default = 600, show_default=True, help = "Number of read pairs to simulate, in millions")
@click.option('-l', '--molecule-length', type = int, default = 100, show_default=True, help = "Mean molecule length (kbp)")
@click.option('-p', '--partitions', type = int, default=1500, show_default=True, help = "How many partitions to generate (×1000)")
@click.option('-m', '--molecules-per', type = int, default = 100, show_default=True, help = "Average number of molecules per partition")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/linkedreads", help = 'Name of output directory')
@click.argument('genome_hap1', required=True, type=click.Path(exists=True), nargs=1)
@click.argument('genome_hap2', required=True, type=click.Path(exists=True), nargs=1)
def reads(genome_hap1, genome_hap2, output_dir, outer_distance, insert_sd, barcodes, read_pairs, molecule_length, partitions, molecules_per, threads, snakemake, quiet, print_only):
    """
    Create linked reads from a genome
 
   
    """


