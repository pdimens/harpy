from .helperfunctions import generate_conda_deps, fetch_file 
from .printfunctions import print_onstart
import rich_click as click
import subprocess
import os
import sys

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/simulate")
@click.option('-f', '--long-frag', type = int, default = 15, show_default= True, help = "Coverage of long fragments")
@click.option('-l', '--long-frag-len', type = int, default = 100, help = "Average length of long fragments (kbp)")
@click.option('-r', '--coverage', type = int, default = 30, help = "Coverage of short reads")
@click.option('-f', '--short-read-len', type = int, default = 150, help = "Length of short reads")
@click.option('-i', '--short-read-insert', type = int, default = 400, help = "Average insert size of short reads")
@click.option('-d', '--insert-sd', type = int, default=50, help = "Standard deviation of insert size (bp)")
@click.option('-m', '--molecules', type = int, default = 16, help = "Average molecules per droplet")
@click.option('-e', '--error-rate', type = float, default = 0.01, help = "Sequencing error rate")
@click.option('-b', '--barcodes', type = click.Path(exists=True), required=True, help = "File of molecular barcodes")
@click.option('-p', '--ploidy', type = click.Choice([1,2]), show_choices = True, default = 2, help = "Ploidy of resulting reads")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/linkedreads", help = 'Name of output directory')
@click.argument('genomes', required=True, type=click.Path(exists=True), nargs=2)
def reads(genomes, output_dir, long_frag, long_frag_len, coverage, short_read_len, short_read_insert, insert_sd, molecules, error_rate, barcodes, ploidy, threads, snakemake, quiet, print_only):
    """
    Create linked reads from a genome
 
   
    """
    # pull sequence error rate files into workflow dir

    # write the parameters into a format-correct parameter file
    with open(f"{workflowdir}/config.txt") as outf:
        outf.write("# haplotypes 1 and 2 of a template genome\n")
        outf.write(f"Path_Fastahap1={genomes[0]}\n")
        outf.write(f"Path_Fastahap2={genomes[1]}\n")
        outf.write("# number of threads to use\n")
        outf.write(f"processors={threads}\n")
        outf.write("# coverage of long fragments\n")
        outf.write(f"CF={coverage}\n")
        outf.write("# average length of long fragments (kbp)\n")
        outf.write(f"Mu_F={long_frag_len}\n")
        outf.write("# length of short reads\n")
        outf.write(f"SR={short_read_len}\n")
        outf.write("# coverage of short reads\n")
        outf.write(f"CR={coverage}\n")
        outf.write("# average insert size of short reads\n")
        outf.write(f"Mu_IS={short_read_insert}\n")
        outf.write("# standard deviation of insert size\n")
        outf.write(f"Std_IS={insert_sd}\n")
        outf.write("# number of molecules per droplet\n")
        outf.write(f"N_FP={molecules}\n")
        outf.write("Fast_mode=N\n")
        outf.write("Seq_error=Y\n")
        outf.write(f"Error_rate={error_rate}\n")
        outf.write("# sequencing error profile\n")
        outf.write(f"Path_Seq_qual=\n")
        outf.write("# barcode error profile\n")
        outf.write(f"Path_Barcode_qual=\n")
        outf.write("# barcode list file\n")
        outf.write(f"Path_barcodepool={barcodes}\n")
        outf.write("# Haploid (1) or Diploid (2)\n")
        outf.write(f"Hap={ploidy}\n")