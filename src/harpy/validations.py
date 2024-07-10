"""Harpy helper functions for validating inputs files to workflows"""

import sys
import os
import re
import gzip
import subprocess
from pathlib import Path
from rich import box, print
from rich.table import Table
import rich_click as click
#from .fileparsers import getnames
from .printfunctions import print_error, print_notice, print_solution, print_solution_with_culprits

def check_envdir(dirpath):
    """Check that the provided dir exists and contains the necessary environment definitions"""
    if not os.path.exists(dirpath):
        print_error("This working directory does not contain the expected directory of conda environment definitions ([blue bold].harpy_envs/[/blue bold])\n  - use [green bold]--conda[/green bold] to recreate it")
        sys.exit(1)
    envlist = os.listdir(dirpath)
    envs = ["align",  "phase", "qc", "r", "simulations", "snp", "stitch", "sv"]
    errcount = 0
    errtable = Table(show_footer=True, box=box.SIMPLE)
    errtable.add_column("File", justify="left", style="blue", no_wrap=True)
    errtable.add_column("Exists", justify="center")
    for i in envs:
        if f"{i}.yaml" in envlist:
            errtable.add_row(f"{i}.yaml", "[blue]✓")
        else:
            errcount += 1
            errtable.add_row(f"{i}.yaml", "[yellow]🗙")
    if errcount > 0:
        print_error(f"The conda environment definition directory ([blue bold]{dirpath}[/blue bold]) is missing [yellow bold]{errcount}[/yellow bold] of the expected definition files. All of the environment files are expected to be present, even if a particular workflow doesn't use it.")
        print_solution_with_culprits(
            "Check that the names conform to Harpy's expectations, otheriwse you can recreate this directory using the [green bold]--conda[/green bold] option.",
            "Expected environment files:"
            )
        print(errtable, file = sys.stderr)
        sys.exit(1)

def validate_input_by_ext(inputfile, option, ext):
    """Check that the input file for a given option has read permissions and matches the acceptable extensions """
    if isinstance(ext, list):
        test = [not(inputfile.lower().endswith(i.lower())) for i in ext]
        if all(test):
            ext_text = " | ".join(ext)
            print_error(f"The input file for [bold]{option}[/bold] must end in one of:\n[green bold]{ext_text}[/green bold]")
            sys.exit(1)
    else:
        if not inputfile.lower().endswith(ext.lower()):
            print_error(f"The input file for [bold]{option}[/bold] must end in [green bold]{ext}[/green bold]")
            sys.exit(1)

def check_impute_params(parameters):
    """Validate the STITCH parameter file for column names, order, types, missing values, etc."""
    with open(parameters, "r", encoding="utf-8") as fp:
        header = fp.readline().rstrip().lower()
        headersplt = header.split()
        correct_header = sorted(["model", "usebx", "bxlimit", "k", "s", "ngen"])
        row = 1
        badrows = []
        badlens = []
        if sorted(headersplt) != correct_header:
            culprits = [i for i in headersplt if i not in correct_header]
            print_error(f"Parameter file [bold]{parameters}[/bold] has incorrect column names. Valid names are:\n[green bold]model usebx bxlimit k s ngen[/green bold]")
            print_solution_with_culprits(
                f"Fix the headers in [bold]{parameters}[/bold] or use [blue bold]harpy stitchparams[/blue bold] to generate a valid parameter file and modify it with appropriate values.",
                "Column names causing this error:"
            )
            click.echo(" ".join(culprits), file = sys.stderr)
            sys.exit(1)
        # instantiate dict with colnames
        data = {}
        for i in headersplt:
            data[i] = []
        while True:
            # Get next line from file
            line = fp.readline()
            row += 1
            # if line is empty, end of file is reached
            if not line:
                break
            if line == "\n":
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
                f"See the problematic rows below. Check that you are using a whitespace (space or tab) delimeter in [bold]{parameters}[/bold] or use [blue bold]harpy stitchparams[/blue bold] to generate a valid parameter file and modify it with appropriate values.",
                "Rows causing this error and their column count:"
            )
            for i in zip(badrows, badlens):
                click.echo(f"{i[0]}\t{i[1]}", file = sys.stderr)
            sys.exit(1)
        
        # Validate each column
        culprits = {}
        colerr = 0
        errtable = Table(show_footer=True, box=box.SIMPLE)
        errtable.add_column("Column", justify="right", style="white", no_wrap=True)
        errtable.add_column("Expected Values", style="blue")
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
            errtable.add_row("usebx", "True, False", ",".join(culprits["usebx"]))
        
        for param in ["bxlimit","k","s","ngen"]:
            culprits[param] = []
            for i,j in enumerate(data[param]):
                if not j.isdigit():
                    culprits[param].append(str(i + 1))
                    colerr += 1
            if culprits[param]:
                errtable.add_row(param, "Integers", ",".join(culprits[param]))

        if colerr > 0:
            print_error(f"Parameter file [bold]{parameters}[/bold] is formatted incorrectly. Not all columns have valid values.")
            print_solution_with_culprits(
                "Review the table below of what values are expected for each column and which rows are causing issues.",
                "Formatting Errors:"
            )
            print(errtable, file = sys.stderr)
            sys.exit(1)

def validate_bam_RG(bamlist):
    """Validate BAM files bamlist to make sure the sample name inferred from the file matches the @RG tag within the file"""
    culpritfiles = []
    culpritIDs   = []
    for i in bamlist:
        samplename = Path(i).stem
        samview = subprocess.run(f"samtools samples {i}".split(), stdout = subprocess.PIPE).stdout.decode('utf-8').split()
        if samplename != samview[0]:
            culpritfiles.append(os.path.basename(i))
            culpritIDs.append(samview[0])
        
    if len(culpritfiles) > 0:
        print_error(f"There are [bold]{len(culpritfiles)}[/bold] alignment files whose ID tags do not match their filenames.")
        print_solution_with_culprits(
            "For alignment files (sam/bam), the base of the file name must be identical to the [green bold]@RD ID:[/green bold] tag in the file header. For example, a file named \'sample_001.bam\' should have the [green bold]@RG ID:sample_001[/green bold] tag. Use the [blue bold]renamebam[/blue bold] script to properly rename alignment files so as to also update the @RG header.",
            "File causing error and their ID tags:"
        )
        for i,j in zip(culpritIDs,culpritfiles):
            click.echo(f"{i}\t{j}", file = sys.stderr)
        sys.exit(1)

def check_phase_vcf(infile):
    """Check to see if the input VCf file is phased or not, govered by the presence of ID=PS or ID=HP tags"""
    vcfheader = subprocess.run(f"bcftools view -h {infile}".split(), stdout = subprocess.PIPE, check = False).stdout.decode('utf-8')
    if ("##FORMAT=<ID=PS" in vcfheader) or ("##FORMAT=<ID=HP" in vcfheader):
        return
    else:
        bn = os.path.basename(infile)
        print_error("The input variant file needs to be phased into haplotypes, but no FORMAT/PS or FORMAT/HP fields were found.")
        print_solution(f"Phase [bold]{bn}[/bold] into haplotypes using [blue bold]harpy phase[/blue bold] or another manner of your choosing and use the phased vcf file as input. If you are confident this file is phased, then the phasing does not follow standard convention and you will need to make sure the phasing information appears as either [bold]FORMAT/PS[/bold] or [bold]FORMAT/HP[/bold] tags.")
        sys.exit(1)

def validate_popfile(infile):
    """Validate the input population file to make sure there are two entries per row"""
    with open(infile, "r", encoding="utf-8") as f:
        rows = [i for i in f.readlines() if i != "\n" and not i.lstrip().startswith("#")]
        invalids = [(i,j) for i,j in enumerate(rows) if len(j.split()) < 2]
        if invalids:
            print_error(f"There are [bold]{len(invalids)}[/bold] rows in [bold]{infile}[/bold] without a space/tab delimiter or don't have two entries for sample[dim]<tab>[/dim]population. Terminating Harpy to avoid downstream errors.")
            print_solution_with_culprits(
                f"Make sure every entry in [bold]{infile}[/bold] uses space or tab delimeters and has both a sample name and population designation. You may comment out rows with a [green bold]#[/green bold] to have Harpy ignore them.",
                "The rows and values causing this error are:"
                )
            _ = [click.echo(f"{i[0]+1}\t{i[1]}", file = sys.stderr) for i in invalids]
            sys.exit(1)

def vcf_samplematch(vcf, bamlist, vcf_samples):
    """Validate that the input VCF file and the samples in the list of BAM files. The directionality of this check is determined by 'vcf_samples', which prioritizes the sample list in the vcf file, rather than bamlist."""
    bcfquery = subprocess.Popen(["bcftools", "query", "-l", vcf], stdout=subprocess.PIPE)
    vcfsamples = bcfquery.stdout.read().decode().split()
    filesamples = [Path(i).stem for i in bamlist]
    if vcf_samples:
        fromthis, query = vcf, vcfsamples
        inthis, search = "the input files", filesamples
    else:
        fromthis, query = "the input files", filesamples
        inthis, search = vcf, vcfsamples

    missing_samples = [x for x in query if x not in search]
    # check that samples in VCF match input directory
    if len(missing_samples) > 0:
        print_error(f"There are [bold]{len(missing_samples)}[/bold] samples found in [bold]{fromthis}[/bold] that are not in [bold]{inthis}[/bold]. Terminating Harpy to avoid downstream errors.")
        print_solution_with_culprits(
            f"[bold]{fromthis}[/bold] cannot contain samples that are absent in [bold]{inthis}[/bold]. Check the spelling or remove those samples from [bold]{fromthis}[/bold] or remake the vcf file to include/omit these samples. Alternatively, toggle [green]--vcf-samples[/green] to aggregate the sample list from the input files or [bold]{vcf}[/bold].",
            "The samples causing this error are:"
        )
        click.echo(", ".join(sorted(missing_samples)) + "\n", file = sys.stderr)
        sys.exit(1)
    return(query)

def validate_popsamples(infiles, popfile, quiet):
    """Validate the presence of samples listed in 'populations' to be in the input files"""
    with open(popfile, "r", encoding="utf-8") as f:
        popsamples = [i.split()[0] for i in f.readlines() if i != "\n" and not i.lstrip().startswith("#")]
    in_samples = [Path(i).stem for i in infiles]
    missing_samples = [x for x in popsamples if x not in in_samples]
    overlooked = [x for x in in_samples if x not in popsamples]
    if len(overlooked) > 0 and not quiet:
        print_notice(f"There are [bold]{len(overlooked)}[/bold] samples found in the inputs that weren\'t included in [bold]{popfile}[/bold]. This will [bold]not[/bold] cause errors and can be ignored if it was deliberate. Commenting or removing these lines will avoid this message. The samples are:\n" + ", ".join(overlooked))
    if len(missing_samples) > 0:
        print_error(f"There are [bold]{len(missing_samples)}[/bold] samples included in [bold]{popfile}[/bold] that weren\'t found in in the inputs. Terminating Harpy to avoid downstream errors.")
        print_solution_with_culprits(
            f"Make sure the spelling of these samples is identical in the inputs and [bold]{popfile}[/bold], or remove them from [bold]{popfile}[/bold].",
            "The samples causing this error are:"
        )
        click.echo(", ".join(sorted(missing_samples)), file = sys.stderr)
        sys.exit(1)

def validate_demuxschema(infile):
    """Validate the file format of the demultiplex schema"""
    with open(infile, "r", encoding="utf-8") as f:
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
        print_solution("Make sure the input file ends with a standard FASTQ extension. These are not case-sensitive.\nAccepted extensions: [green bold].fq .fastq .fq.gz .fastq.gz[/green bold]")
        sys.exit(1)
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
        sys.exit(1)

def validate_regions(regioninput, genome):
    """validates the --regions input of harpy snp to infer whether it's an integer, region, or file"""
    try:
        # is an int
        region = int(regioninput)
    except:
        region = regioninput
    if isinstance(region, int) and region < 10:
        print_error("Input for [green bold]--regions[/green bold] was interpreted as an integer to create equal windows to call variants. Integer input for [green bold]--regions[/green bold] must be greater than or equal to [blue bold]10[/blue bold].")
        sys.exit(1)
    elif isinstance(region, int) and region >= 10:
        return "windows"
    else:
        pass
    # is a string
    reg = re.split(r"[\:-]", regioninput)
    if len(reg) == 3:
        # is a single region, check the types to be [str, int, int]
        err = ""
        try:
            reg[1] = int(reg[1])
        except:
            err += f"The region start position [blue bold]({reg[1]})[/blue bold] is not a valid integer. "
        try:
            reg[2] = int(reg[2])
        except:
            err += f"The region end position [blue bold]({reg[2]})[/blue bold] is not a valid integer."
        if err != "":
            print_error("Input for [green bold]--regions[/green bold] was interpreted as a single region. " + err)
            print_solution("If providing a single region to call variants, it should be in the format [yellow bold]contig:start-end[/yellow bold], where [yellow bold]start[/yellow bold] and [yellow bold]end[/yellow bold] are integers. If the input is a file and was incorrectly interpreted as a region, try changing the name to avoid using colons ([yellow bold]:[/yellow bold]) or dashes ([yellow bold]-[/yellow bold]).")
            sys.exit(1)
        # check if the region is in the genome

        contigs = {}
        if genome.lower().endswith("gz"):
            with gzip.open(genome, "r") as fopen:
                for line in fopen:
                    line = line.decode()
                    if line.startswith(">"):
                        cn = line.rstrip("\n").lstrip(">").split()[0]
                        contigs[cn] = 0
                    else:
                        contigs[cn] += len(line.rstrip("\n")) - 1
        else:
            with open(genome, "r", encoding="utf-8") as fout:
                for line in fopen:
                    if line.startswith(">"):
                        cn = line.rstrip("\n").lstrip(">").split()[0]
                        contigs[cn] = 0
                    else:
                        contigs[cn] += len(line.rstrip("\n")) - 1
        err = ""
        if reg[0] not in contigs:
            print_error(f"The contig ([bold yellow]{reg[0]})[/bold yellow]) of the input region [yellow bold]{regioninput}[/yellow bold] was not found in [bold]{genome}[/bold].")
            sys.exit(1)
        if reg[1] > contigs[reg[0]]:
            err += f"- start position: [bold yellow]{reg[1]}[/bold yellow]\n"
        if reg[2] > contigs[reg[0]]:
            err += f"- end position: [bold yellow]{reg[2]}[/bold yellow]"
        if err != "":
            print_error(f"Components of the input region [yellow bold]{regioninput}[/yellow bold] were not found in [bold]{genome}[/bold]:\n" + err)
            sys.exit(1)
        return "region"
    # is a file specifying regions
    if not os.path.isfile(regioninput):
        print_error("Input for [green bold]--regions[/green bold] was interpreted as a file of regions to call variants over and this file was not found.")
        print_solution(f"Check that the path to [bold]{regioninput}[/bold] is correct.")
        sys.exit(1)
    with open(regioninput, "r", encoding="utf-8") as fin:
        badrows = []
        idx = 0
        while True:
            line = fin.readline()
            if not line:
                break
            else:
                idx += 1
            row = line.split()
            if len(row) != 3:
                badrows.append(idx)
            else:
                try:
                    int(row[1])
                    int(row[2])
                except:
                    badrows.append(idx)

        if badrows:
            print_error("The input file is formatted incorrectly.")
            print_solution_with_culprits(
                "Rows in the provided file need to be [bold]space[/bold] or [bold]tab[/bold] delimited with the format [yellow bold]contig start end[/yellow bold] where [yellow bold]start[/yellow bold] and [yellow bold]end[/yellow bold] are integers.",
                "Rows triggering this error:"
                )
            click.echo(",".join([i for i in badrows]), file = sys.stderr)
            sys.exit(1)
    return "file"