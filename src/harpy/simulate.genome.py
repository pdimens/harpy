from .helperfunctions import generate_conda_deps, fetch_file 
from .printfunctions import print_onstart
import rich_click as click
import subprocess
import os
import sys

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/simulate")
@click.option('-t', '--variant-type', type = click.Choice(["snp, indel", "inversion", "cnv", "translocation"], help = "Type of variant to simulate"))
@click.option('-v', '--vcf', type=click.Path(exists=True), help = 'VCF file of known variants to simulate')
@click.option('-h', '--heterozygosity', type = click.FloatRange(0,1), default = 0, help = "Percent heterozygosity if simulating diploid")
@click.option('-p', '--parameters', type = click.Path(exists=True), help = "Simulation parameter file")
@click.option('-e', '--exclude-chr', type = click.Path(exists=True), help = "Text file of chromosomes to ignore")
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/genome", show_default=True, help = 'Name of output directory')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=1)
def genome(input, output_dir, snps, indels, inversions, cnv, translocations, parameters, exclude_chr, snakemake, quiet, print_only):
    """
    Introduce variants into a haploid genome
 
    You can create a haploid genome with variants introduced into it or diploid genome
    by setting a heterozygosity value greater than `0`. Use either a VCF file to simulate
    known variants or provide a YAML file for `--parameters` to simulate random variants
    of that type. The parameter YAML file can be created using `harpy simparams` and
    modified accordingly. The simulator (simuG) can only simulate one type of variant
    at a time and you will need to run this module subsequent times to simulate additional
    kinds of variants.
    """
