from .helperfunctions import generate_conda_deps, fetch_file 
#from .fileparsers import get_samples_from_fastq, parse_fastq_inputs
from .printfunctions import print_onstart
import rich_click as click
import subprocess
import os
import sys

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/qc")
@click.option('-s', '--snps', type=click.Path(exists=True) help = 'VCF file of SNPs to simulate')
@click.option('-i', '--indels', type=click.Path(exists=True), help = 'VCF file of indels to simulate')
@click.option('-i', '--inversions', type = type=click.Path(exists=True), help = 'VCF file of inversions to simulate')
@click.option('-c', '--cnv', type=click.Path(exists=True), help = 'VCF file of CNVs to simulate')
@click.option('-t', '--translocations', type=click.Path(exists=True), help = 'VCF file of translocations to simulate')
@click.option('-h', '--heterozygosity', type = float, help = "Percent heterozygosity between both haplotypes")
@click.option('-p', '--parameters', type = click.Path(exists=True), help = "Simulation parameter file")
@click.option('-e', '--exclude-chr', type = click.Path(exists=True), help = "Text file of chromosomes to ignore")
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/genome", show_default=True, help = 'Name of output directory')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def genome(input, output_dir, snps, indels, inversions, cnv, translocations, parameters, exclude_chr, snakemake, quiet, print_only):
    """
    Introduce variants into a haploid genome
 
    You can create a haploid genome with variants introduced into it or diploid genome
    by setting a heterozygosity value. By using any one of the command line arguments
     requiring a VCF file, you will introduce known variants of that type
    into the provided input genome. However, if you want to randomly simulate variants, that
    requires more granular configuration by way of the `--parameters` file, which can be created
    using `harpy simparams` and modified accordingly.
    """
