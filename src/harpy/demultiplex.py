import rich_click as click
import subprocess
import sys
import os

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

#TODO demux sample sheet validation

@click.command(no_args_is_help = True)
@click.option('-f', '--file', required = True, type=click.Path(exists=True), metavar = "File Path", help = 'The forward (or reverse) multiplexed FASTQ file')
@click.option('-b', '--samplesheet', required = True, type=click.Path(exists=True), metavar = "File Path", help = 'Tab-delimited file of sample<tab>barcode')
@click.option('-m', '--method', default = "gen1", show_default = True, type = click.Choice(["gen1"], case_sensitive = False), metavar = "String", help = "Haplotag technology of the sequences")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def demultiplex(file, method, samplesheet, threads, snakemake, quiet):
    """
    Demultiplex haplotagged FASTQ files

    Use one of the four gzipped FASTQ files provided by the sequencer (I1, I2, R1, R2).fastq.gz for the `--file` argument, Harpy will infer the other three.
    Double-check that you are using the correct haplotag method (`--method`), since the different barcoding approaches
    have very different demultiplexing strategies. Note: the `--samplesheet` must be or space delimited and have no header (i.e. no column names).
    """
    with open(samplesheet, "r") as f:
        #rows = [i.split("\t")[0] for i in f.readlines()]
        badrows = []
        for i,j in enumerate(f.readlines()):
            # a casual way to ignore empty lines or lines with >2 fields
            try:
                smpl, bc = j.split()
                d[smpl] = bc
            except:
                badrows.append((i,j))
                continue
    if badrows:
        #TODO ERROR MESSAGE HERE JUST LIKE POPCHECK

    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/demultiplex.{method}.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"infile={file}")
    command.append(f"samplefile={samplesheet}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)