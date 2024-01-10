from .helperfunctions import getnames, createregions, validate_bamfiles, validate_popfile, validate_vcfsamples
import rich_click as click
import subprocess
import sys
import os

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-g', '--genome', type=click.Path(exists=True), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-p', '--populations', type=click.Path(exists = True), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, metavar = "Integer", help = 'Ploidy of samples')
@click.option('-w', '--windowsize', default = 50000, show_default = True, type = int, metavar = "Integer", help = "Interval size for parallel variant calling")
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def mpileup(genome, threads, directory, populations, ploidy, windowsize, extra_params, snakemake, quiet, print_only):
    """
    Call variants from using bcftools mpileup
    
    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`. Available methods are:
    """
    samplenames = getnames(directory, '.bam')
    callcoords, linkedgenome = createregions(genome, windowsize, "mpileup")
    directory = directory.rstrip("/^")
    validate_bamfiles(directory, samplenames)
    command = (f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/variants-mpileup.smk').split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"seq_directory={directory}")
    command.append(f"samplenames={samplenames}")
    if populations is not None:
        rows = validate_popfile(populations)
        # check that samplenames and populations line up
        validate_vcfsamples(directory, populations, samplenames, rows, quiet)
        command.append(f"groupings={populations}")
    command.append(f"ploidy={ploidy}")
    command.append(f"windowsize={windowsize}")
    command.append(f"intervals={callcoords}")
    command.append(f"genomefile={linkedgenome}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    if print_only:
        click.echo(" ".join(command))
    else:
        _module = subprocess.run(command)
        sys.exit(_module.returncode)

@click.command(no_args_is_help = True)
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-x', '--ploidy', default = 2, show_default = True, type=int, metavar = "Integer", help = 'Ploidy of samples')
@click.option('-w', '--windowsize', default = 50000, show_default = True, type = int, metavar = "Integer", help = "Interval size for parallel variant calling")
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def freebayes(genome, threads, directory, populations, ploidy, windowsize, extra_params, snakemake, quiet, print_only):
    """
    Call variants using freebayes
    
    Optionally specify `--populations` for population-aware variant calling.
    Use **harpy popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    samplenames = getnames(directory, '.bam')
    callcoords, linkedgenome = createregions(genome, windowsize, "freebayes")
    directory = directory.rstrip("/^")
    validate_bamfiles(directory, samplenames)
    command = (f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/variants-freebayes.smk').split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"seq_directory={directory}")
    command.append(f"samplenames={samplenames}")
    if populations is not None:
        rows = validate_popfile(populations)
        # check that samplenames and populations line up
        validate_vcfsamples(directory, populations, samplenames, rows, quiet)
        command.append(f"groupings={populations}")
    command.append(f"ploidy={ploidy}")
    command.append(f"windowsize={windowsize}")
    command.append(f"intervals={callcoords}")
    command.append(f"genomefile={linkedgenome}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    if print_only:
        click.echo(" ".join(command))
    else:
        _module = subprocess.run(command)
        sys.exit(_module.returncode)