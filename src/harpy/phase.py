from .helperfunctions import fetch_file, generate_conda_deps, getnames, vcfcheck, vcf_samplematch, validate_bamfiles
import sys
import os
import subprocess
import rich_click as click

@click.command(no_args_is_help = True)
@click.option('-v', '--vcf', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Path to BCF/VCF file')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance threshold to separate molecules')
@click.option('-p', '--prune-threshold', default = 7, show_default = True, type = click.IntRange(0,100), metavar = "Integer", help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-b', '--ignore-bx',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Ignore barcodes when phasing')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Use samples present in vcf file for phasing rather than those found the directory')
@click.option('-g', '--genome', type=click.Path(exists=True, dir_okay=False), metavar = "File path", help = 'Path to genome assembly if wanting to also extract reads spanning indels')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional HapCut2 parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 2, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake',  type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def phase(vcf, directory, threads, molecule_distance, prune_threshold, vcf_samples, genome, snakemake, extra_params, ignore_bx, quiet, print_only):
    """
    Phase SNPs into haplotypes

    You may choose to omit barcode information with `--ignore-bx`, although it's usually
    better to include that information. Use `--vcf-samples` to phase only
    the samples present in your input `--vcf` file rather than all the samples present in
    the `--directory`.
    """
    fetch_file("phase-pop.smk", "Phase/workflow/")
    fetch_file("HapCut2.Rmd", "Phase/workflow/report/")
    directory = directory.rstrip("/^")
    vcfcheck(vcf)
   # if vcf.lower().endswith(".vcf.gz"):
   #     click.echo(f"Notice: HapCut2 does not accept gzipped vcf files. Converting to bcf.")
   #     variantfile = vcf[0:-7] + ".bcf"
   #     subprocess.run(f"bcftools view {vcf} -Ob > {variantfile}".split())
   # else:
   #     variantfile = vcf
    
    samplenames = vcf_samplematch(vcf, directory, vcf_samples)
    validate_bamfiles(directory, samplenames)
    prune_threshold /= 100
    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake --cores {threads} --directory . --snakefile Phase/workflow/phase-pop.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    if genome is not None:
        command.append(f"indels={genome}")
        if not os.path.exists(f"{genome}.fai"):
            subprocess.run(f"samtools faidx --fai-idx {genome}.fai {genome}".split())
    command.append(f"seq_directory={directory}")
    command.append(f"samplenames={samplenames}")
    command.append(f"variantfile={vcf}")
    command.append(f"noBX={ignore_bx}")
    command.append(f"prune={prune_threshold}")
    command.append(f"molecule_distance={molecule_distance}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    if print_only:
        click.echo(" ".join(command))
    else:
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)