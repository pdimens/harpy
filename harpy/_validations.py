"""Harpy helper functions for validating inputs files to workflows"""

import sys
import os
import re
import gzip
import subprocess
from pathlib import Path
from rich import box
from rich import print as rprint
from rich.table import Table
import rich_click as click
from ._printing import print_error, print_notice, print_solution, print_solution_with_culprits
from ._misc import harpy_progressbar
from concurrent.futures import ThreadPoolExecutor, as_completed

def is_gzip(file_path):
    """helper function to determine if a file is gzipped"""
    try:
        with gzip.open(file_path, 'rt') as f:
            f.read(10)
        return True
    except gzip.BadGzipFile:
        return False

def is_plaintext(file_path):
    """helper function to determine if a file is plaintext"""
    try:
        with open(file_path, 'r') as f:
            f.read(10)
        return True
    except UnicodeDecodeError:
        return False

def check_envdir(dirpath):
    """Check that the provided dir exists and contains the necessary environment definitions"""
    if not os.path.exists(dirpath):
        print_error("missing conda files", "This working directory does not contain the expected directory of conda environment definitions ([blue bold].harpy_envs/[/blue bold])\n  - use [green bold]--conda[/green bold] to recreate it")
        sys.exit(1)
    envlist = os.listdir(dirpath)
    envs = ["align", "metassembly", "phase", "qc", "r", "simulations", "stitch", "variants"]
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
        print_error("missing conda files", f"The conda environment definition directory ([blue bold]{dirpath}[/blue bold]) is missing [yellow bold]{errcount}[/yellow bold] of the expected definition files. All of the environment files are expected to be present, even if a particular workflow doesn't use it.")
        print_solution_with_culprits(
            "Check that the names conform to Harpy's expectations, otheriwse you can recreate this directory using the [green bold]--conda[/green bold] option.",
            "Expected environment files:"
            )
        rprint(errtable, file = sys.stderr)
        sys.exit(1)

def validate_input_by_ext(inputfile, option, ext):
    """Check that the input file for a given option has read permissions and matches the acceptable extensions """
    if isinstance(ext, list):
        test = [not(inputfile.lower().endswith(i.lower())) for i in ext]
        if all(test):
            ext_text = " | ".join(ext)
            print_error("invalid file extension", f"The input file for [bold]{option}[/bold] must end in one of:\n[green bold]{ext_text}[/green bold]")
            sys.exit(1)
    else:
        if not inputfile.lower().endswith(ext.lower()):
            print_error("invalid file extension", f"The input file for [bold]{option}[/bold] must end in [green bold]{ext}[/green bold]")
            sys.exit(1)

def check_impute_params(parameters):
    """Validate the STITCH parameter file for column names, order, types, missing values, etc."""
    with open(parameters, "r", encoding="utf-8") as paramfile:
        header = paramfile.readline().rstrip().lower()
        headersplt = header.split()
        colnames = ["name", "model", "usebx", "bxlimit", "k", "s", "ngen"]
        correct_header = sorted(colnames)
        n_cols = len(colnames)
        row = 1
        badrows = []
        badlens = []
        if sorted(headersplt) != correct_header:
            invalid_col = [i for i in headersplt if i not in colnames]
            missing_col = [i for i in colnames if i not in headersplt]
            errtext = []
            culprit_text = []
            if invalid_col:
                errtext.append("\n  - invalid column names")
                culprit_text.append(f"[red]Invalid columns:[/red] {", ".join(invalid_col)}")
            if missing_col:
                errtext.append("\n  - missing columns")
                culprit_text.append( f"[yellow]Missing columns:[/yellow] {", ".join(missing_col)}")

            print_error("incorrect columns", f"Parameter file [bold]{parameters}[/bold] has incorrect column names{"".join(errtext)}\nValid names are: [green bold]{" ".join(colnames)}[/green bold]")
            print_solution_with_culprits(
                f"Fix the headers in [bold]{parameters}[/bold] or use [blue bold]harpy stitchparams[/blue bold] to generate a valid parameter file and modify it with appropriate values.",
                "Column causing this error:"
            )
            rprint("\n".join(culprit_text), file = sys.stderr)
            sys.exit(1)
        else:
            # in case the columns are out of order, reorder the columns to match `colnames`
            col_order = [colnames.index(item) for item in headersplt]
        # instantiate dict with colnames
        data = {}
        while True:
            # Get next line from file
            line = paramfile.readline()
            row += 1
            # if line is empty, end of file is reached
            if not line:
                break
            if line == "\n":
                # be lenient with empty rows
                continue
            # split the line by whitespace and reorder to match expected colname order
            row_values = [line.rstrip().split()[i] for i in col_order]
            if len(row_values) == n_cols: 
                data[row_values[0]] = dict(zip(colnames[1:], row_values[1:]))
            else:
                badrows.append(row)
                badlens.append(len(row_values))
        if len(badrows) > 0:
            print_error("invalid rows", f"Parameter file [blue]{parameters}[/blue] is formatted incorrectly. Not all rows have the expected {n_cols} columns.")
            print_solution_with_culprits(
                f"See the problematic rows below. Check that you are using a whitespace (space or tab) delimeter in [blue]{parameters}[/blue] or use [blue green]harpy stitchparams[/blue green] to generate a valid parameter file and modify it with appropriate values.",
                "Rows causing this error and their column count:"
            )
            for i in zip(badrows, badlens):
                click.echo(f"{i[0]}\t{i[1]}", file = sys.stderr)
            sys.exit(1)
        
        # validate each row
        row_error = False
        errtable = Table(show_footer=True, box=box.SIMPLE)
        errtable.add_column("Row", justify="right", style="white", no_wrap=True)
        errtable.add_column("Columns with Issues", style = "white")

        row = 1
        for k,v in data.items():
            badcols = []
            if re.search(r'[^a-z0-9\-_\.]',k, flags = re.IGNORECASE):
                badcols.append("name")
            if v["model"] not in ["pseudoHaploid", "diploid","diploid-inbred"]:
                badcols.append("model")
            if f"{v["usebx"]}".lower() not in ["true", "false", "yes", "y", "no", "n"]:
                badcols.append("usebx")
            else:
                if f"{v["usebx"]}".lower() in ["true", "yes", "y"]:                
                    v["usebx"] = True
                else:
                    v["usebx"] = False
            for param in ["bxlimit", "k", "s", "ngen"]:
                if not v[param].isdigit():
                    badcols.append(param)
                else:
                    v[param] = int(v[param])
            if badcols:
                row_error = True
                errtable.add_row(f"{row}", ", ".join(badcols))
            row += 1
        if row_error:
            print_error("invalid parameter values", f"Parameter file [bold]{parameters}[/bold] is formatted incorrectly. Some rows have incorrect values for one or more parameters.")
            print_solution_with_culprits(
                "Review the table below of which rows/columns are causing issues",
                "Formatting Errors:"
            )
            rprint(errtable, file = sys.stderr)
            sys.exit(1)        
        return data

def validate_bam_RG(bamlist, threads, quiet):
    """Validate BAM files bamlist to make sure the sample name inferred from the file matches the @RG tag within the file"""
    culpritfiles = []
    culpritIDs   = []
    def check_RG(bamfile):
        samplename = Path(bamfile).stem
        samview = subprocess.run(f"samtools samples {bamfile}".split(), stdout = subprocess.PIPE).stdout.decode('utf-8').split()
        if samplename != samview[0]:
            return os.path.basename(i), samview[0]

    with harpy_progressbar(quiet) as progress:
        task_progress = progress.add_task("[green]Checking RG tags...", total=len(bamlist))
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(check_RG, i) for i in bamlist]
            for future in as_completed(futures):
                result = future.result()
                progress.update(task_progress, advance=1)
                if result:
                    culpritfiles.append(result[0])
                    culpritIDs.append(result[1])
        
    if len(culpritfiles) > 0:
        print_error("sample ID mismatch",f"There are [bold]{len(culpritfiles)}[/bold] alignment files whose RG tags do not match their filenames.")
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
        print_error("vcf not phased", "The input variant file needs to be phased into haplotypes, but no [green]FORMAT/PS[/green] or [green]FORMAT/HP[/green] fields were found.")
        print_solution(f"Phase [bold]{bn}[/bold] into haplotypes using [blue bold]harpy phase[/blue bold] or another manner of your choosing and use the phased vcf file as input. If you are confident this file is phased, then the phasing does not follow standard convention and you will need to make sure the phasing information appears as either [green]FORMAT/PS[/green] or [green]FORMAT/HP[/green] tags.")
        sys.exit(1)

def validate_popfile(infile):
    """Validate the input population file to make sure there are two entries per row"""
    with open(infile, "r", encoding="utf-8") as f:
        rows = [i for i in f.readlines() if i != "\n" and not i.lstrip().startswith("#")]
        invalids = [(i,j) for i,j in enumerate(rows) if len(j.split()) < 2]
        if invalids:
            print_error("invalid format", f"There are [bold]{len(invalids)}[/bold] rows in [bold]{infile}[/bold] without a space/tab delimiter or don't have two entries for sample[dim]<tab>[/dim]population. Terminating Harpy to avoid downstream errors.")
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
        print_error("mismatched inputs", f"There are [bold]{len(missing_samples)}[/bold] samples found in [blue]{fromthis}[/blue] that are not in [blue]{inthis}[/blue]. Terminating Harpy to avoid downstream errors.")
        print_solution_with_culprits(
            f"[blue]{fromthis}[/blue] cannot contain samples that are absent in [blue]{inthis}[/blue]. Check the spelling or remove those samples from [blue]{fromthis}[/blue] or remake the vcf file to include/omit these samples. Alternatively, toggle [green]--vcf-samples[/green] to aggregate the sample list from the input files or [blue]{vcf}[/blue].",
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
        print_notice(f"There are [bold]{len(overlooked)}[/bold] samples found in the inputs that weren\'t included in [blue]{popfile}[/blue]. This will [bold]not[/bold] cause errors and can be ignored if it was deliberate. Commenting or removing these lines will avoid this message. The samples are:\n" + ", ".join(overlooked))
    if len(missing_samples) > 0:
        print_error("mismatched inputs", f"There are [bold]{len(missing_samples)}[/bold] samples included in [blue]{popfile}[/blue] that weren\'t found in in the inputs. Terminating Harpy to avoid downstream errors.")
        print_solution_with_culprits(
            f"Make sure the spelling of these samples is identical in the inputs and [blue]{popfile}[/blue], or remove them from [blue]{popfile}[/blue].",
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
            print_error(f"invalid format", "There are [bold]{len(invalids)}[/bold] rows in [blue]{infile}[/blue] without a space/tab delimiter or don't have two entries for sample[dim]<tab>[/dim]barcode. Terminating Harpy to avoid downstream errors.")
            print_solution_with_culprits(
                f"Make sure every entry in [blue]{infile}[/blue] uses space or tab delimeters and has both a sample name and barcode designation. You may comment out rows with a [green]#[/green] to have Harpy ignore them.",
                "The rows and values causing this error are:"
                )
            _ = [click.echo(f"{i[0]+1}\t{i[1]}", file = sys.stderr) for i in invalids]
            sys.exit(1)

def check_demux_fastq(file):
    """Check for the presence of corresponding FASTQ files from a single provided FASTQ file based on pipeline expectations."""
    bn = os.path.basename(file)
    if not bn.lower().endswith("fq") and not bn.lower().endswith("fastq") and not bn.lower().endswith("fastq.gz") and not bn.lower().endswith("fq.gz"):     
        print_error("unrecognized extension", f"The file [blue]{bn}[/blue] is not recognized as a FASTQ file by the file extension.")
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
        print_error("missing files", f"Not all necessary files with prefix [bold]{prefix}[/bold] present")
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
        print_error("window size too small", "Input for [green bold]--regions[/green bold] was interpreted as an integer to create equal windows to call variants. Integer input for [green bold]--regions[/green bold] must be greater than or equal to [green bold]10[/green bold].")
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
            err += f"The region start position [green bold]({reg[1]})[/green bold] is not a valid integer. "
        try:
            reg[2] = int(reg[2])
        except:
            err += f"The region end position [green bold]({reg[2]})[/green bold] is not a valid integer."
        if err != "":
            print_error("invalid region", "Input for [green bold]--regions[/green bold] was interpreted as a single region. " + err)
            print_solution("If providing a single region to call variants, it should be in the format [yellow bold]contig:start-end[/yellow bold], where [yellow bold]start[/yellow bold] and [yellow bold]end[/yellow bold] are integers. If the input is a file and was incorrectly interpreted as a region, try changing the name to avoid using colons ([yellow bold]:[/yellow bold]) or dashes ([yellow bold]-[/yellow bold]).")
            sys.exit(1)
        # check if the region is in the genome

        contigs = {}
        if is_gzip(genome):
            with gzip.open(genome, "rt") as fopen:
                for line in fopen:
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
            print_error("contig not found", f"The contig ([bold yellow]{reg[0]})[/bold yellow]) of the input region [yellow bold]{regioninput}[/yellow bold] was not found in [blue]{genome}[/blue].")
            sys.exit(1)
        if reg[1] > contigs[reg[0]]:
            err += f"- start position: [bold yellow]{reg[1]}[/bold yellow]\n"
        if reg[2] > contigs[reg[0]]:
            err += f"- end position: [bold yellow]{reg[2]}[/bold yellow]"
        if err != "":
            print_error("region out of bounds", f"Components of the input region [yellow bold]{regioninput}[/yellow bold] were not found in [blue]{genome}[/blue]:\n" + err)
            sys.exit(1)
        return "region"
    # is a file specifying regions
    if not os.path.isfile(regioninput):
        print_error("file not found", "Input for [green bold]--regions[/green bold] was interpreted as a file and this file was not found.")
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
            print_error("invalid format", "The input file is formatted incorrectly.")
            print_solution_with_culprits(
                "Rows in the provided file need to be [bold]space[/bold] or [bold]tab[/bold] delimited with the format [yellow bold]contig start end[/yellow bold] where [yellow bold]start[/yellow bold] and [yellow bold]end[/yellow bold] are integers.",
                "Rows triggering this error:"
                )
            click.echo(",".join([i for i in badrows]), file = sys.stderr)
            sys.exit(1)
    return "file"

def check_fasta(genofile, quiet):
    """perform validations on fasta file for extensions and file contents"""
    ext_options = [".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn"]
    ext_correct = 0
    for i in ext_options:
        if genofile.lower().endswith(i) or genofile.lower().endswith(i + ".gz"):
            ext_correct += 1
    if ext_correct == 0 and not quiet:
        print_notice(f"[blue]{genofile}[/blue] has an unfamiliar FASTA file extension. Common FASTA file extensions are: [green]" + ", ".join(ext_options) + "[/green] and may also be gzipped.")

    # validate fasta file contents
    if is_gzip(genofile):
        fasta = gzip.open(genofile, 'rt')
    elif is_plaintext(genofile):
        fasta = open(genofile, 'r', encoding="utf-8")
    else:
        print_error("unknown file type", f"Unable to determine file encoding for [blue]{genofile}[/blue]. Please check that it is a gzipped or uncompressed FASTA file.")
        sys.exit(1)
    line_num = 0
    seq_id = 0
    seq = 0
    last_header = False
    for line in fasta:
        line_num += 1
        if line.startswith(">"):
            seq_id += 1
            if last_header:
                print_error("consecutive contig names", f"All contig names must be followed by at least one line of nucleotide sequences, but two consecutive lines of contig names were detected. This issue was identified at line [bold]{line_num}[/bold] in [blue]{genofile}[/blue], but there may be others further in the file.")
                print_solution("See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
                sys.exit(1)
            else:
                last_header = True
            if len(line.rstrip()) == 1:
                print_error("unnamed contigs", f"All contigs must have an alphanumeric name, but a contig was detected without a name. This issue was identified at line [bold]{line_num}[/bold] in [blue]{genofile}[/blue], but there may be others further in the file.")
                print_solution("See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
                sys.exit(1)
            if line.startswith("> "):
                print_error("invalid contig names", f"All contig names must be named [green bold]>contig_name[/green bold], without a space, but a contig was detected with a space between the [green bold]>[/green bold] and contig_name. This issue was identified at line [bold]{line_num}[/bold] in [blue]{genofile}[/blue], but there may be others further in the file.")
                print_solution("See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
                sys.exit(1)
        elif line == "\n":
            print_error("empty lines", f"Empty lines are not permitted in FASTA files, but one was detected at line [bold]{line_num}[/bold] in [blue]{genofile}[/blue]. The scan ended at this error, but there may be others further in the file.")
            print_solution("See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
            sys.exit(1)
        else:
            seq += 1
            last_header = False
    fasta.close()
    solutiontext = "FASTA files must have at least one contig name followed by sequence data on the next line. Example:\n"
    solutiontext += "[green]  >contig_name\n  ATACAGGAGATTAGGCA[/green]\n"
    # make sure there is at least one of each
    if seq_id == 0:
        print_error("contig names absent", f"No contig names detected in [blue]{genofile}[/blue].")
        print_solution(solutiontext + "\nSee the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
        sys.exit(1)
    if seq == 0:
        print_error("sequences absent", f"No sequences detected in [blue]{genofile}[/blue].")
        print_solution(solutiontext + "\nSee the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
        sys.exit(1)


def validate_fastq_bx(fastq_list, threads, quiet):
    def validate(fastq):
        BX = False
        BC = False
        if is_gzip(fastq):
            fq = gzip.open(fastq, "rt")
        else:
            fq = open(fastq, "r")
        with fq:
            while True:
                line = fq.readline()
                if not line:
                    break
                if not line.startswith("@"):
                    continue
                BX = True if "BX:Z" in line else BX
                BC = True if "BC:Z" in line else BC
                if BX and BC:
                    print_error("clashing barcode tags", f"Both [green bold]BC:Z[/green bold] and [green bold]BX:Z[/green bold] tags were detected in the read headers for [blue]{os.path.basename(fastq)}[/blue]. Athena accepts [bold]only[/bold] one of [green bold]BC:Z[/green bold] or [green bold]BX:Z[/green bold].")
                    print_solution("Check why your data has both tags in use and remove/rename one of the tags.")
                    sys.exit(1)
            # check for one or the other after parsing is done
            if not BX and not BC:
                print_error("no barcodes found",f"No [green bold]BC:Z[/green bold] or [green bold]BX:Z[/green bold] tags were detected in read headers for [blue]{os.path.basename(fastq)}[/blue]. Athena requires the linked-read barcode to be present as either [green bold]BC:Z[/green bold] or [/green bold]BX:Z[/green bold] tags.")
                print_solution("Check that this is linked-read data and that the barcode is demultiplexed from the sequence line into the read header as either a `BX:Z` or `BC:Z` tag.")
                sys.exit(1)

    # parellelize over the fastq list
    with harpy_progressbar(quiet) as progress:
        task_progress = progress.add_task("[green]Validating FASTQ inputs...", total=len(fastq_list))
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(validate, i) for i in fastq_list]
            for future in as_completed(futures):
                progress.update(task_progress, advance=1)