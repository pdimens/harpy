from .harpymisc import getnames_err, vcfcheck, check_impute_params
import rich_click as click
import subprocess
import sys
import os

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

#@click.option('-f', '--filter', is_flag=True, help="Filter VCF file to keep SNPs with QUAL>20 and DP>10")
@click.command(no_args_is_help = True)
@click.option('-v', '--vcf', required = True, type=click.Path(exists=True),metavar = "File Path", help = 'Path to BCF/VCF file')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-p', '--parameters', required = True, type=click.Path(exists=True), metavar = "File Path", help = 'STITCH parameter file (tab-delimited)')
#@click.option('-x', '--extra-params', default = "", type = str, metavar = "String", help = 'Additional STITCH parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Use samples present in vcf file for imputation rather than those found the directory')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def impute(parameters, directory, threads, vcf, vcf_samples, snakemake, quiet):
    """
    Impute genotypes using variants and sequences
    
    Requires a parameter file, use **harpy extra --stitch-params** to generate one and adjust it for your study.
    """
    ## validate inputs ##
    vcfcheck(vcf)
    samplenames = getnames_err(directory, '.bam')
    ### check that samples in VCF match input directory
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

    check_impute_params(parameters)

    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/impute.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    # generate and store list of viable contigs (minimum of 2 biallelic SNPs)
    # doing it here so it doesn't have to run each time inside the workflow
    vbn = os.path.basename(vcf)
    if not os.path.exists(f"Impute/input/_{vbn}.list"):
        os.makedirs("Impute/input/", exist_ok = True)
        print("\033[1mPreprocessing:\033[00m Identifying contigs with at least 2 biallelic SNPs", file = sys.stderr)
        biallelic = subprocess.Popen(f"bcftools view -M2 -v snps {vcf} -Ob".split(), stdout = subprocess.PIPE)
        contigs = subprocess.Popen("""bcftools query -f '%CHROM\\n'""".split(), stdin = biallelic.stdout, stdout = subprocess.PIPE)
        c_sort = subprocess.Popen("sort", stdin = contigs.stdout, stdout = subprocess.PIPE)
        unq = subprocess.Popen("uniq -c".split(), stdin = c_sort.stdout, stdout = subprocess.PIPE)
        contigs_out = subprocess.run(["awk", r'{ if ($1 > 1) {print $2} }'], stdin = unq.stdout, stdout = subprocess.PIPE).stdout.splitlines()
        contigs = []
        with open(f"Impute/input/_{vbn}.list", "w") as f:
            for l in contigs_out:
                l_corr = l.decode().replace("\'", "")
                _ = f.write(f"{l_corr}\n")
                contigs.append(f"{l_corr}")
    else:
        with open(f"Impute/input/_{vbn}.list", "r") as f:
            contigs = [line.rstrip() for line in f]
    command.append('--config')
    directory = directory.rstrip("/^")
    command.append(f"seq_directory={directory}")
    command.append(f"samplenames={samplenames}")
    command.append(f"variantfile={vcf}")
    command.append(f"paramfile={parameters}")
    command.append(f"contigs={contigs}")
    #command.append(f"extra={extra_params}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)