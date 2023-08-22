from .harpymisc import getnames_err, vcfcheck
import rich_click as click
import subprocess
import sys
import os

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-v', '--vcf', required = True, type=click.Path(exists=True), metavar = "File Path", help = 'Path to BCF/VCF file')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-m', '--molecule-distance', default = 50000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance delineating separate molecules')
@click.option('-p', '--prune-threshold', default = 7, show_default = True, type = click.IntRange(0,100), metavar = "Integer", help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-b', '--ignore-bx',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Ignore barcodes when phasing')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Use samples present in vcf file for phasing rather than those found the directory')
@click.option('-g', '--genome', type=click.Path(exists=True), metavar = "File path", help = 'Path to genome assembly if wanting to also extract reads spanning indels')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional HapCut2 parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 2, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake',  type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def phase(vcf, directory, threads, molecule_distance, prune_threshold, vcf_samples, genome, snakemake, extra_params, ignore_bx, quiet):
    """
    Phase SNPs into haplotypes

    You may choose to omit barcode information with `--ignore-bx`, although it's usually
    better to include that information. Use the `--vcf-samples` toggle to phase only
    the samples present in your input `--vcf` file rather than all the samples present in
    the `--directory`.
    """
    vcfcheck(vcf)
    if vcf.lower().endswith(".vcf.gz"):
        print(f"Notice: HapCut2 does not accept gzipped vcf files. Converting to bcf.")
        variantfile = vcf[0:-7] + ".bcf"
        subprocess.run(f"bcftools view {vcf} -Ob > {variantfile}".split())
    else:
        variantfile = vcf
    
    bcfquery = subprocess.Popen(["bcftools", "query", "-l", vcf], stdout=subprocess.PIPE)
    if vcf_samples:
        samplenames = bcfquery.stdout.read().decode().split()
        s_list = getnames_err(directory, '.bam')
        fromthis = vcf
        inthis = directory
    else:
        samplenames = getnames_err(directory, '.bam')
        s_list = bcfquery.stdout.read().decode().split()
        fromthis = directory
        inthis = vcf
    missing_samples = [x for x in samplenames if x not in s_list]
    # check that samples in VCF match input directory
    if len(missing_samples) > 0:
        print(f"\n\033[1;33mERROR:\033[00m There are {len(missing_samples)} samples found in \033[01m{fromthis}\033[00m that are not in \033[01m{inthis}\033[00m. Terminating Harpy to avoid downstream errors. The samples causing this error are:", file = sys.stderr)
        print(", ".join(sorted(missing_samples)), file = sys.stderr)
        print(f"\n\033[1;34mSOLUTION:\033[00m \033[01m{fromthis}\033[00m cannot contain samples that are absent in \033[01m{inthis}\033[00m. Check the spelling or remove those samples from \033[01m{fromthis}\033[00m or remake the vcf file to include/omit these samples. Alternatively, toggle \033[01m--vcf-samples\033[00m to aggregate the sample list from \033[01m{directory}\033[00m or \033[01m{vcf}\033[00m.\n", file = sys.stderr)
        sys.exit(1)
    prune_threshold /= 100
    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/phase-pop.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    directory = directory.rstrip("/^")
    if indels is not None:
        command.append(f"indels={indels}")
        if not os.path.exists(indels + ".fai"):
            sys.run(f"samtools faidx {indels}".split())
    command.append(f"seq_directory={directory}")
    command.append(f"samplenames={samplenames}")
    command.append(f"variantfile={variantfile}")
    command.append(f"noBX={ignore_bx}")
    command.append(f"prune={prune_threshold}")
    command.append(f"molecule_distance={molecule_distance}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)