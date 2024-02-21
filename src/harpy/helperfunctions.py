import sys
import os
import re
import shutil
import glob
import gzip
import subprocess
from collections import Counter
from urllib.request import urlretrieve
from pathlib import Path
from rich import print
from rich.panel import Panel
from rich.table import Table
import rich_click as click

## define some rich print functions for less redundancy
def print_onstart(text, title):
    """Print a panel of info on workflow run"""
    click.echo("")
    print(Panel(text, title = f"[bold]Harpy {title}", title_align = "left", border_style = "white", subtitle= "Initializing"), file = sys.stderr)
def print_error(errortext):
    """Print a yellow panel with error text"""
    print(Panel(errortext, title = "[bold]Error", title_align = "left", border_style = "yellow"), file = sys.stderr)

def print_solution(solutiontext):
    """Print a blue panel with solution text"""
    print(Panel(solutiontext, title = "[bold]Solution", title_align = "left", border_style = "blue"), file = sys.stderr)

def print_solution_with_culprits(solutiontext, culprittext):
    """Print a blue panel with solution text and the list of offenders below it"""
    print(Panel(solutiontext, title = "[bold]Solution", title_align = "left", subtitle = culprittext, border_style = "blue"), file = sys.stderr)

def print_notice(noticetext):
    """Print a white panel with information text text"""
    print(Panel(noticetext, title = "Notice", title_align = "left", border_style = "white"), file = sys.stderr)

## recurring checks and such ##
def vcfcheck(vcf):
    """Check that a file ends with one of .vcf, .bcf, or .vcf.gz. Is case insensitive."""
    vfile = vcf.lower()
    if vfile.endswith(".vcf") or vfile.endswith(".bcf") or vfile.endswith(".vcf.gz"):
        pass
    else:
        print_error(f"Supplied variant call file [bold]{vcf}[/bold] must end in one of [.vcf | .vcf.gz | .bcf]")
        exit(1)

def getnames(directory, ext):
    """Find all files in 'directory' that end with 'ext'"""
    samplenames = set([i.split(ext)[0] for i in os.listdir(directory) if i.endswith(ext)])
    if len(samplenames) < 1:
        print_error(f"No sample files ending with [bold]{ext}[/bold] found in [bold]{directory}[/bold].")
        sys.exit(1)
    return samplenames

def get_samples_from_fastq(directory):
    """Identify the sample names from a directory containing FASTQ files"""
    full_flist = [i for i in glob.iglob(f"{directory}/*") if not os.path.isdir(i)]
    r = re.compile(r".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    full_fqlist = list(filter(r.match, full_flist))
    fqlist = [os.path.basename(i) for i in full_fqlist]
    bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
    if len(fqlist) == 0:
        print_error(f"No fastq files with acceptable names found in [bold]{directory}[/bold]")
        print_solution("Check that the file endings conform to [green].[/green][[green]F[/green][dim]|[/dim][green]R1[/green]][green].[/green][[green]fastq[/green][dim]|[/dim][green]fq[/green]][green].gz[/green]\nRead the documentation for details: https://pdimens.github.io/harpy/haplotagdata/#naming-conventions")
        sys.exit(1)

    return set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

def createregions(infile, window, method):
    """Create a BED file of genomic intervals of size 'window'. Uses 1- or 0- based numbering depending on mpileup or freebayes 'method'"""
    bn = os.path.basename(infile)
    os.makedirs("Genome", exist_ok = True)
    base = 0 if method == "freebayes" else 1
    gen_zip = True if bn.lower().endswith(".gz") else False
    if method == "freebayes":
        # freebayes requires uncompressed genome
        if gen_zip:
            # remove .gz extension
            bn = bn[:-3]
            if not os.path.exists(f"Genome/{bn}"):
                with open(f"Genome/{bn}", "w") as fo:
                    subprocess.run(f"gzip -dc {infile}".split(), stdout = fo)
        else:
            if not os.path.exists(f"Genome/{bn}"):
                subprocess.run(f"ln -sr {infile} Genome/{bn}".split())
    else:
        if not os.path.exists(f"Genome/{bn}"):
            ftype = subprocess.run(["file", infile], stdout=subprocess.PIPE).stdout.decode('utf-8')
            if "Blocked GNU Zip" in ftype:
                # is bgzipped, just link it
                subprocess.run(f"ln -sr {infile} Genome/{bn}".split())
            elif "gzip compressed data" in ftype:
                # is regular gzipped, needs to be bgzipped
                subprocess.run(f"zcat {infile} | bgzip -c > Genome/{bn}".split())
            else:
                # not compressed, just link
                subprocess.run(f"ln -sr {infile} Genome/{bn}".split())

    if not os.path.exists(f"Genome/{bn}.fai"):
        try:
            subprocess.run(f"samtools faidx --fai-idx Genome/{bn}.fai --gzi-idx Genome/{bn}.gzi Genome/{bn}".split(), stderr = subprocess.DEVNULL)
        except:
            subprocess.run(f"samtools faidx --fai-idx Genome/{bn}.fai Genome/{bn}".split(), stderr = subprocess.DEVNULL)

    with open(f"Genome/{bn}.fai") as fai:
        bedregion = []
        while True:
            # Get next line from file
            line = fai.readline()
            # if line is empty, end of file is reached
            if not line:
                break
            # split the line by tabs
            lsplit = line.split()
            contig = lsplit[0]
            c_len = int(lsplit[1])
            c_len = c_len - 1 if base == 0 else c_len
            start = base
            end = window
            starts = [base]
            ends = [window]
            while end < c_len:
                end = end + window if (end + window) < c_len else c_len
                ends.append(end)
                start += window
                starts.append(start)
            for (startpos, endpos) in zip (starts,ends):
                bedregion.append(f"{contig}:{startpos}-{endpos}")
        return bedregion, f"Genome/{bn}"

def check_impute_params(parameters):
    """Validate the STITCH parameter file for column names, order, types, missing values, etc."""
    with open(parameters, "r") as fp:
        header = fp.readline().rstrip().lower()
        headersplt = header.split()
        correct_header = sorted(["model", "usebx", "bxlimit", "k", "s", "ngen"])
        row = 1
        badrows = []
        badlens = []
        if sorted(headersplt) != correct_header:
            culprits = [i for i in headersplt if i not in correct_header]
            print_error(f"Parameter file [bold]{parameters}[/bold] has incorrect column names. Valid names are:\n[green]model usebx bxlimit k s ngen[/green]")
            print_solution_with_culprits(
                f"Fix the headers in [bold]{parameters}[/bold] or use [green]harpy stitchparams[/green] to generate a valid parameter file and modify it with appropriate values.",
                "Column names causing this error:"
            )
            click.echo(" ".join(culprits), file = sys.stderr)
            sys.exit(1)
        # instantiate dict with colnames
        data = dict()
        for i in headersplt:
            data[i] = []
        while True:
            # Get next line from file
            line = fp.readline()
            row += 1
            # if line is empty, end of file is reached
            if not line:
                break
            # split the line by whitespace
            rowvals = line.rstrip().split()
            rowlen = len(rowvals)
            if rowlen != 6:
                badrows.append(row)
                badlens.append(rowlen)
            elif len(badrows) > 0:
                continue
            else:
                for j,k in zip(headersplt, rowvals):
                    data[j].append(k)

        if len(badrows) > 0:
            print_error(f"Parameter file [bold]{parameters}[/bold] is formatted incorrectly. Not all rows have the expected 6 columns.")
            print_solution_with_culprits(
                f"See the problematic rows below. Check that you are using a whitespace (space or tab) delimeter in [bold]{parameters}[/bold] or use [green]harpy stitchparams[/green] to generate a valid parameter file and modify it with appropriate values.",
                "Rows causing this error and their column count:"
            )
            for i in zip(badrows, badlens):
                click.echo(f"{i[0]}\t{i[1]}", file = sys.stderr)
            sys.exit(1)
        
        # Validate each column
        culprits = dict()
        colerr = 0
        errtable = Table(title="Formatting Errors")
        errtable.add_column("Column", justify="right", style="white", no_wrap=True)
        errtable.add_column("Expected Values", style="green")
        errtable.add_column("Rows with Issues", style = "white")

        culprits["model"] = []
        for i,j in enumerate(data["model"]):
            if j not in ["pseudoHaploid", "diploid","diploid-inbred"]:
                culprits["model"].append(str(i + 1))
                colerr += 1
        if culprits["model"]:
            errtable.add_row("model", "diploid, diploid-inbred, pseudoHaploid", ", ".join(culprits["model"]))

        culprits["usebx"] = []
        for i,j in enumerate(data["usebx"]):
            if j not in [True, "TRUE", "true", False, "FALSE", "false", "Y","y", "YES", "Yes", "yes", "N", "n", "NO", "No", "no"]:
                culprits["usebx"].append(str(i + 1))
                colerr += 1
        if culprits["usebx"]:
            errtable.add_row("usebx", "True, False", ", ".join(culprits["usebx"]))
        
        for param in ["bxlimit","k","s","ngen"]:
            culprits[param] = []
            for i,j in enumerate(data[param]):
                if not j.isdigit():
                    culprits[param].append(str(i + 1))
                    colerr += 1
            if culprits[param]:
                errtable.add_row(param, "Integers", ", ".join(culprits[param]))

        if colerr > 0:
            print_error(f"Parameter file [bold]{parameters}[/bold] is formatted incorrectly. Not all columns have valid values.")
            print_solution("Review the table below of what values are expected for each column and which rows are causing issues.")
            print(errtable, file = sys.stderr)
            exit(1)

def validate_bamfiles(dir, namelist):
    """Validate BAM files in directory 'dir' to make sure the sample name inferred from the file matches the @RG tag within the file"""
    culpritfiles = []
    culpritIDs   = []
    for i in namelist:
        fname = f"{dir}/{i}.bam"
        samview = subprocess.Popen(f"samtools view -H {fname}".split(), stdout = subprocess.PIPE)
        IDtag = subprocess.run("grep ^@RG".split(), stdin = samview.stdout, stdout=subprocess.PIPE).stdout.decode('utf-8')
        r = re.compile("(\t)(ID:.*?)(\t)")
        IDsearch = r.search(IDtag)
        if IDsearch:
            # strip right and left whitespaces and rm "ID:"
            IDval = IDsearch.group(0).rstrip().lstrip()[3:]
            # does the ID: tag match the sample name?
            if IDval != i:
                culpritfiles.append(os.path.basename(fname))
                culpritIDs.append(IDval)
        else:
            culpritfiles.append(os.path.basename)
            culpritIDs.append("missing @RG ID:")
        
    if len(culpritfiles) > 0:
        print_error(f"There are [bold]{len(culpritfiles)}[/bold] alignment files whose ID tags do not match their filenames.")
        print_solution_with_culprits(
            f"For alignment files (sam/bam), the base of the file name must be identical to the [green]@RD ID:[/green] tag in the file header. For example, a file named \'sample_001.bam\' should have the [green]@RG ID:sample_001[/green] tag. Use the [green]renamebam[/green] script to properly rename alignment files so as to also update the @RG header.",
            "File causing error and their ID tags:"
        )
        for i,j in zip(culpritfiles,culpritIDs):
            click.echo(f"{i}\t{j}", file = sys.stderr)
        exit(1)

def check_phase_vcf(infile):
    """Check to see if the input VCf file is phased or not, govered by the presence of ID=PS or ID=HP tags"""
    vcfheader = subprocess.run(f"bcftools view -h {infile}".split(), stdout = subprocess.PIPE).stdout.decode('utf-8')
    if ("##FORMAT=<ID=PS" in vcfheader) or ("##FORMAT=<ID=HP" in vcfheader):
        return
    else:
        bn = os.path.basename(infile)
        print_error(f"The input variant file needs to be phased into haplotypes, but no FORMAT/PS or FORMAT/HP fields were found.")
        print_solution(f"Phase [bold]{bn}[/bold] into haplotypes using [green]harpy phase[/green] or another manner of your choosing and use the phased vcf file as input. If you are confident this file is phased, then the phasing does not follow standard convention and you will need to make sure the phasing information appears as either [bold]FORMAT/PS[/bold] or [bold]FORMAT/HP[/bold] tags.")
        exit(1)

def validate_popfile(infile):
    """Validate the input population file to make sure there are two entries per row"""
    with open(infile, "r") as f:
        rows = [i for i in f.readlines() if i != "\n" and not i.lstrip().startswith("#")]
        invalids = [(i,j) for i,j in enumerate(rows) if len(j.split()) < 2]
        if invalids:
            print_error(f"There are [bold]{len(invalids)}[/bold] rows in [bold]{infile}[/bold] without a space/tab delimiter or don't have two entries for sample[dim]<tab>[/dim]population. Terminating Harpy to avoid downstream errors.")
            print_solution_with_culprits(
                f"Make sure every entry in [bold]{infile}[/bold] uses space or tab delimeters and has both a sample name and population designation. You may comment out rows with a [green]#[/green] to have Harpy ignore them.",
                "The rows and values causing this error are:"
                )
            _ = [click.echo(f"{i[0]+1}\t{i[1]}", file = sys.stderr) for i in invalids]
            sys.exit(1)
        else:
            return rows

def vcf_samplematch(vcf, directory, vcf_samples):
    """Validate that the input VCF file and the samples in the 'directory' (BAM files). The directionality of this check is determined by 'vcf_samples', which prioritizes the sample list in the file, rather that directory."""
    bcfquery = subprocess.Popen(["bcftools", "query", "-l", vcf], stdout=subprocess.PIPE)
    filesamples = bcfquery.stdout.read().decode().split()
    dirsamples  = getnames(directory, '.bam')
    if vcf_samples:
        fromthis, query = vcf, filesamples
        inthis, search = directory, dirsamples
    else:
        fromthis, query = directory, dirsamples
        inthis, search = vcf, filesamples

    missing_samples = [x for x in query if x not in search]
    # check that samples in VCF match input directory
    if len(missing_samples) > 0:
        print_error(f"There are [bold]{len(missing_samples)}[/bold] samples found in [bold]{fromthis}[/bold] that are not in [bold]{inthis}[/bold]. Terminating Harpy to avoid downstream errors.")
        print_solution_with_culprits(
            f"[bold]{fromthis}[/bold] cannot contain samples that are absent in [bold]{inthis}[/bold]. Check the spelling or remove those samples from [bold]{fromthis}[/bold] or remake the vcf file to include/omit these samples. Alternatively, toggle [green]--vcf-samples[/green] to aggregate the sample list from [bold]{directory}[/bold] or [bold]{vcf}[/bold].",
            "The samples causing this error are:"
        )
        click.echo(", ".join(sorted(missing_samples)), file = sys.stderr)
        sys.exit(1)
    return(query)

def validate_vcfsamples(directory, populations, samplenames, rows, quiet):
    """Validate the presence of samples listed in 'populations' to be in the target directory"""
    p_list = [i.split()[0] for i in rows]
    missing_samples = [x for x in p_list if x not in samplenames]
    overlooked = [x for x in samplenames if x not in p_list]
    if len(overlooked) > 0 and not quiet:
        print_notice(f"There are [bold]{len(overlooked)}[/bold] samples found in [bold]{directory}[/bold] that weren\'t included in [bold]{populations}[/bold]. This will not cause errors and can be ignored if it was deliberate. Commenting or removing these lines will avoid this message. The samples are:\n" + ", ".join(overlooked))
    if len(missing_samples) > 0:
        print_error(f"There are [bold]{len(missing_samples)}[/bold] samples included in [bold]{populations}[/bold] that weren\'t found in [bold]{directory}[/bold]. Terminating Harpy to avoid downstream errors.")
        print_solution_with_culprits(
            f"Make sure the spelling of these samples is identical in [bold]{directory}[/bold] and [bold]{populations}[/bold], or remove them from [bold]{populations}[/bold].",
            "The samples causing this error are:"
        )
        click.echo(", ".join(sorted(missing_samples)), file = sys.stderr)
        sys.exit(1)

def validate_demuxschema(infile):
    """Validate the file format of the demultiplex schema"""
    with open(infile, "r") as f:
        rows = [i for i in f.readlines() if i != "\n" and not i.lstrip().startswith("#")]
        invalids = [(i,j) for i,j in enumerate(rows) if len(j.split()) < 2]
        if invalids:
            print_error(f"There are [bold]{len(invalids)}[/bold] rows in [bold]{infile}[/bold] without a space/tab delimiter or don't have two entries for sample[dim]<tab>[/dim]barcode. Terminating Harpy to avoid downstream errors.")
            print_solution_with_culprits(
                f"Make sure every entry in [bold]{infile}[/bold] uses space or tab delimeters and has both a sample name and barcode designation. You may comment out rows with a [green]#[/green] to have Harpy ignore them.",
                "The rows and values causing this error are:"
                )
            _ = [click.echo(f"{i[0]+1}\t{i[1]}", file = sys.stderr) for i in invalids]
            sys.exit(1)

def check_demux_fastq(file):
    """Check for the presence of corresponding FASTQ files from a single provided FASTQ file based on pipeline expectations."""
    bn = os.path.basename(file)
    if not bn.lower().endswith("fq") and not bn.lower().endswith("fastq") and not bn.lower().endswith("fastq.gz") and not bn.lower().endswith("fq.gz"):     
        print_error(f"The file {bn} is not recognized as a FASTQ file by the file extension.")
        print_solution("Make sure the input file ends with a standard FASTQ extension. These are not case-sensitive.\nAccepted extensions: [blue].fq .fastq .fq.gz .fastq.gz[/blue]")
        exit(1)
    ext = re.search(r"(?:\_00[0-9])*\.f(.*?)q(?:\.gz)?$", file, re.IGNORECASE).group(0)
    prefix     = re.sub(r"[\_\.][IR][12]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", "", bn)
    prefixfull = re.sub(r"[\_\.][IR][12]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", "", file)
    filelist = []
    printerr = False
    for i in ["I1", "I2","R1","R2"]:
        chkfile = f"{prefixfull}_{i}{ext}"
        TF = os.path.exists(chkfile)
        printerr = True if not TF else printerr
        symbol = " " if TF else "X"
        filelist.append(f"\033[91m{symbol}\033[0m  {prefix}_{i}{ext}")
    if printerr:
        print_error(f"Not all necessary files with prefix [bold]{prefix}[/bold] present")
        print_solution_with_culprits(
            "Demultiplexing requires 4 sequence files: the forward and reverse sequences (R1 and R2), along with the forward and reverse indices (I2 and I2).",
            "Necessary/expected files:"
        )
        _ = [click.echo(i, file = sys.stderr) for i in filelist]
        exit(1)

def generate_conda_deps():
    """Create the YAML files of the workflow conda dependencies"""
    condachannels = ["conda-forge", "bioconda", "defaults"]
    environ = {
        "qc" : ["falco", "fastp"],
        "align": ["bwa", "ema","icu","libzlib", "samtools=1.19", "seqtk", "xz"],
        "variants.snp": ["bcftools=1.19", "freebayes=1.3.6"],
        "variants.sv": ["leviathan", "naibr-plus"],
        "phase" : ["hapcut2", "whatshap"],
        "r-env" : ["bioconductor-complexheatmap", "r-circlize", "r-dt", "r-flexdashboard", "r-ggplot2", "r-ggridges", "r-plotly", "r-tidyr", "r-stitch"]
    }

    os.makedirs("harpyenvs", exist_ok = True)

    for i in environ:
        # don't overwrite existing
        if not os.path.isfile(f"harpyenvs/{i}.yaml"):
            with open(f"harpyenvs/{i}.yaml", "w") as yml:
                yml.write(f"name: {i}\n")
                yml.write("channels:\n  - ")
                yml.write("\n  - ".join(condachannels))
                yml.write("\ndependencies:\n  - ")
                yml.write("\n  - ".join(environ[i]) + "\n")


def fetch_file(file, destination, rename=None):
    """Find the 'file' in the PATH and copy it to the 'destination' with the original metadata"""
    result = None
    try:
        result = subprocess.check_output(["whereis", file])
    except:
        print_error(f"The GNU program \'whereis\', which is used to locate the file, was not found on the system and therefore unable to retrieve it. Terminating harpy.")
        print_solution("Make sure \'whereis\' is installed on the system. It is usually provided by default in all Unix-like operating systems.")

    if result is None:
        return []

    result = result.decode().splitlines()
    for line in result:
        if line.endswith(":"):
            print_error(f"The file \"{file}\" was not found in PATH, cannot run Harpy module.")
            print_solution(f"Make sure harpy was installed correctly and that you are in the harpy conda environment.")
            click.echo("See documentation: https://pdimens.github.io/harpy/install/")
            exit(1)
        else:
            result = line.split(" ")[-1]
            break
    
    os.makedirs(destination, exist_ok = True)
    if rename:
        destination += rename
    # copy2 to keep metadata during copy
    shutil.copy2(result, destination)

def biallelic_contigs(vcf):
    """Identify which contigs have at least 2 biallelic SNPs"""
    vbn = os.path.basename(vcf)
    if not os.path.exists(f"Impute/input/_{vbn}.list"):
        os.makedirs("Impute/input/", exist_ok = True)
        biallelic = subprocess.Popen(f"bcftools view -M2 -v snps {vcf} -Ob".split(), stdout = subprocess.PIPE)
        contigs = subprocess.run("""bcftools query -f '%CHROM\\n'""".split(), stdin = biallelic.stdout, stdout = subprocess.PIPE).stdout.decode().splitlines()
        counts = Counter(contigs)
        contigs = [i.replace("\'", "") for i in counts if counts[i] > 1]
        with open(f"Impute/input/_{vbn}.list", "w") as f:
            _ = [f.write(f"{i}\n") for i in contigs]
    else:
        with open(f"Impute/input/_{vbn}.list", "r") as f:
            contigs = [line.rstrip() for line in f]
    
    if len(contigs) == 0:
        print_error("No contigs with at least 2 biallelic SNPs identified. Cannot continue with imputation.")
        exit(1)
    else:
        return contigs