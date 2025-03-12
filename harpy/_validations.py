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
from ._misc import harpy_progressbar, safe_read
from concurrent.futures import ThreadPoolExecutor, as_completed

# logic to properly refresh progress bar for jupyter sessions
try:
    __IPYTHON__
    _in_ipython_session = True
except NameError:
    _in_ipython_session = False

def is_gzip(file_path: str) -> bool:
    """helper function to determine if a file is gzipped"""
    try:
        with gzip.open(file_path, 'rt') as f:
            f.read(10)
        return True
    except gzip.BadGzipFile:
        return False

def is_plaintext(file_path: str) -> bool:
    """helper function to determine if a file is plaintext"""
    try:
        with open(file_path, 'r') as f:
            f.read(10)
        return True
    except UnicodeDecodeError:
        return False

def check_impute_params(parameters: str) -> dict:
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
                culprit_text.append(f"[red]Invalid columns:[/red] " + ", ".join(invalid_col))
            if missing_col:
                errtext.append("\n  - missing columns")
                culprit_text.append( f"[yellow]Missing columns:[/yellow] " + ", ".join(missing_col))

            print_error("incorrect columns", f"Parameter file [bold]{parameters}[/bold] has incorrect column names" + "".join(errtext) + f"\nValid names are: [green bold]" + " ".join(colnames) + "[/green bold]")
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
            if v["usebx"].lower() not in ["true", "false", "yes", "y", "no", "n"]:
                badcols.append("usebx")
            else:
                if v["usebx"].lower() in ["true", "yes", "y"]:                
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

def validate_bam_RG(bamlist: str, threads: int, quiet: int) -> None:
    """Validate BAM files bamlist to make sure the sample name inferred from the file matches the @RG tag within the file"""
    culpritfiles = []
    culpritIDs   = []
    def check_RG(bamfile):
        samplename = Path(bamfile).stem
        samview = subprocess.run(f"samtools samples {bamfile}".split(), stdout = subprocess.PIPE).stdout.decode('utf-8').split()
        if samplename != samview[0]:
            return os.path.basename(bamfile), samview[0]

    with harpy_progressbar(quiet) as progress, ThreadPoolExecutor(max_workers=threads) as executor:
        task_progress = progress.add_task("[green]Checking RG tags", total=len(bamlist))
        futures = [executor.submit(check_RG, i) for i in bamlist]
        for future in as_completed(futures):
            result = future.result()
            progress.advance(task_progress)
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

def check_phase_vcf(infile: str) -> None:
    """Check to see if the input VCf file is phased or not, govered by the presence of ID=PS or ID=HP tags"""
    vcfheader = subprocess.run(f"bcftools view -h {infile}".split(), stdout = subprocess.PIPE, check = False).stdout.decode('utf-8')
    if ("##FORMAT=<ID=PS" in vcfheader) or ("##FORMAT=<ID=HP" in vcfheader):
        return
    else:
        bn = os.path.basename(infile)
        print_error("vcf not phased", "The input variant file needs to be phased into haplotypes, but no [green]FORMAT/PS[/green] or [green]FORMAT/HP[/green] fields were found.")
        print_solution(f"Phase [bold]{bn}[/bold] into haplotypes using [blue bold]harpy phase[/blue bold] or another manner of your choosing and use the phased vcf file as input. If you are confident this file is phased, then the phasing does not follow standard convention and you will need to make sure the phasing information appears as either [green]FORMAT/PS[/green] or [green]FORMAT/HP[/green] tags.")
        sys.exit(1)

def validate_popfile(infile: str) -> None:
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

def vcf_sample_match(vcf: str, bamlist: list[str], vcf_samples: bool) -> list[str]:
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

def vcf_contig_match(contigs: list[str], vcf: str) -> None:
    vcf_contigs = []
    vcf_header = subprocess.run(["bcftools", "head", vcf], capture_output = True, text = True)
    for i in vcf_header.stdout.split("\n"):
        if i.startswith("##contig="):
            unformatted_contig = i.split(",")[0]
            vcf_contigs.append(unformatted_contig.lstrip("##contig=<ID="))
    bad_names = []
    for i in contigs:
        if i not in vcf_contigs:
            bad_names.append(i)
    if bad_names:
        shortname = os.path.basename(vcf)
        print_error("contigs absent", f"Some of the provided contigs were not found in [blue]{shortname}[/blue]. This will definitely cause plotting errors in the workflow.")
        print_solution_with_culprits("Check that your contig names are correct, including uppercase and lowercase. It's possible that you listed a contig in the genome that isn't in the variant call file due to filtering.", f"Contigs absent in {shortname}:")
        click.echo(",".join([i for i in bad_names]), file = sys.stderr)
        sys.exit(1)

def validate_popsamples(infiles: list[str], popfile: str, quiet: int) -> None:
    """Validate the presence of samples listed in 'populations' to be in the input files"""
    with open(popfile, "r", encoding="utf-8") as f:
        popsamples = [i.split()[0] for i in f.readlines() if i != "\n" and not i.lstrip().startswith("#")]
    in_samples = [Path(i).stem for i in infiles]
    missing_samples = [x for x in popsamples if x not in in_samples]
    overlooked = [x for x in in_samples if x not in popsamples]
    if len(overlooked) > 0 and not quiet != 2:
        print_notice(f"There are [bold]{len(overlooked)}[/bold] samples found in the inputs that weren\'t included in [blue]{popfile}[/blue]. This will [bold]not[/bold] cause errors and can be ignored if it was deliberate. Commenting or removing these lines will avoid this message. The samples are:\n" + ", ".join(overlooked))
    if len(missing_samples) > 0:
        print_error("mismatched inputs", f"There are [bold]{len(missing_samples)}[/bold] samples included in [blue]{popfile}[/blue] that weren\'t found in in the inputs. Terminating Harpy to avoid downstream errors.")
        print_solution_with_culprits(
            f"Make sure the spelling of these samples is identical in the inputs and [blue]{popfile}[/blue], or remove them from [blue]{popfile}[/blue].",
            "The samples causing this error are:"
        )
        click.echo(", ".join(sorted(missing_samples)), file = sys.stderr)
        sys.exit(1)

def validate_demuxschema(infile:str, return_len: bool = False) -> None | int:
    """Validate the file format of the demultiplex schema. Set return_len to True to return the number of samples"""
    code_letters = set() #codes can be Axx, Bxx, Cxx, Dxx
    segment_ids = set()
    samples = set()
    segment_pattern = re.compile(r'^[A-D]\d{2}$')
    with open(infile, 'r') as file:
        for line in file:
            try:
                sample, segment_id = line.rstrip().split()
                if not segment_pattern.match(segment_id):
                    print_error("invalid segment format", f"Segment ID [green]{segment_id}[/green] does not follow the expected format.")
                    print_solution("This haplotagging design expects segments to follow the format of letter [green bold]A-D[/green bold] followed by [bold]two digits[/bold], e.g. [green bold]C51[/green bold]). Check that your ID segments or formatted correctly and that you are attempting to demultiplex with a workflow appropriate for your data design.")
                    sys.exit(1)
                code_letters.add(segment_id[0])
                samples.add(sample)
                if segment_id in segment_ids:
                    print_error("ambiguous segment ID", "An ID segment must only be associated with a single sample.")
                    print_solution_with_culprits(
                        "A barcode segment can only be associated with a single sample. For example: [green bold]C05[/green bold] cannot identify both [green]sample_01[/green] and [green]sample_2[/green]. In other words, a segment can only appear once.",
                        "The segment triggering this error is:"
                        )
                    click.echo(segment_id)
                    sys.exit(1)
                else:
                    segment_ids.add(segment_id)
            except ValueError:
                # skip rows without two columns
                continue
    if not code_letters:
        print_error("incorrect schema format", f"Schema file [blue]{os.path.basename(infile)}[/blue] has no valid rows. Rows should be sample<tab>segment, e.g. sample_01<tab>C75")
        sys.exit(1)
    if len(code_letters) > 1:
        print_error("invalid schema", f"Schema file [blue]{os.path.basename(infile)}[/blue] has sample IDs occurring  in different barcode segments.")
        print_solution_with_culprits(
            "All sample IDs for this barcode design should be in a single segment, such as [bold green]C[/bold green] or [bold green]D[/bold green]. Make sure the schema contains only one segment.",
            "The segments identified in the schema:"
        )
        click.echo(", ".join(code_letters))
        sys.exit(1)
    if return_len:
        return len(samples)

def validate_regions(regioninput: int | str, genome: str) -> str:
    """validates the --regions input of harpy snp to infer whether it's an integer, region, or file"""
    try:
        # is an int
        region = int(regioninput)
        if region < 10:
            print_error("window size too small", "Input for [green bold]--regions[/green bold] was interpreted as an integer to create equal windows to call variants. Integer input for [green bold]--regions[/green bold] must be greater than or equal to [green bold]10[/green bold].")
            sys.exit(1)
        else:
            return "windows"
    except ValueError:
        region = regioninput
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
        with safe_read(genome) as fopen:
            for line in fopen:
                if line.startswith(">"):
                    cn = line.rstrip("\n").lstrip(">").split()[0]
                    contigs[cn] = 0
                else:
                    contigs[cn] += len(line.rstrip("\n"))
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
        for idx, line in enumerate(fin, 1):
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

def check_fasta(genofile: str) -> str:
    """perform validations on fasta file for extensions and file contents"""
    # validate fasta file contents
    line_num = 0
    seq_id = 0
    seq = 0
    last_header = False
    with safe_read(genofile) as fasta:
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

def fasta_contig_match(contigs: str, fasta: str) -> None:
    """Checks whether a list of contigs are present in a fasta file"""
    valid_contigs = []
    with safe_read(fasta) as gen_open:
        for line in gen_open:
            if not line.startswith(">"):
                continue
            # get just the name, no comments or starting character
            valid_contigs.append(line.rstrip().lstrip(">").split()[0])
    bad_names = []
    for i in contigs:
        if i not in valid_contigs:
            bad_names.append(i)
    if bad_names:
        shortname = os.path.basename(fasta)
        print_error("contigs absent", f"Some of the provided contigs were not found in [blue]{shortname}[/blue]. This will definitely cause plotting errors in the workflow.")
        print_solution_with_culprits("Check that your contig names are correct, including uppercase and lowercase.", f"Contigs absent in {shortname}:")
        click.echo(",".join([i for i in bad_names]), file = sys.stderr)
        sys.exit(1)

def validate_fastq_bx(fastq_list: list[str], threads: int, quiet: int) -> None:
    """
    Parse a list of fastq files to verify that they have BX/BC tag, and only one of those two types per file
    """
    def validate(fastq):
        BX = False
        BC = False
        with safe_read(fastq) as fq:
            for line in fq:
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
    with harpy_progressbar(quiet) as progress, ThreadPoolExecutor(max_workers=threads) as executor:
        task_progress = progress.add_task("[green]Validating FASTQ inputs", total=len(fastq_list))
        futures = [executor.submit(validate, i) for i in fastq_list]
        for future in as_completed(futures):
            progress.advance(task_progress)

def validate_barcodefile(infile: str, return_len: bool = False, quiet: int = 0, limit: int = 60) -> None | int:
    """Does validations to make sure it's one length, within a length limit, one per line, and nucleotides"""
    if is_gzip(infile):
        print_error("incorrect format", f"The input file must be in uncompressed format. Please decompress [blue]{infile}[/blue] and try again.")
        sys.exit(1)
    lengths = set()
    nucleotides = {'A','C','G','T'}
    def validate(line_num, bc_line):
        barcode = bc_line.rstrip()
        if len(barcode.split()) > 1:
            print_error("incorrect format", f"There must be one barcode per line, but multiple entries were detected on [bold]line {line_num}[/bold] in [blue]{infile}[/blue]")
            sys.exit(1)
        if not set(barcode).issubset(nucleotides) or barcode != barcode.upper():
            print_error("incorrect format", f"Invalid barcode format on [bold]line {line_num }[/bold]: [yellow]{barcode}[/yellow].\nBarcodes in [blue]{infile}[/blue] must be captial letters and only contain standard nucleotide characters [green]ATCG[/green].")
            sys.exit(1)
        return len(barcode)

    with open(infile, "r") as bc_file, harpy_progressbar(quiet) as progress:
        out = subprocess.Popen(['wc', '-l', infile], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
        linenum = int(out.partition(b' ')[0])
        if linenum > 96**4:
            print_error("Too many barcodes", f"The maximum number of barcodes possible with haplotagging is [bold]96^4 (84,934,656)[/bold], but there are [yellow]{linenum}[/yellow] barcodes in [blue]{infile}[/blue]. Please use fewer barcodes.")
            sys.exit(1)
        task_progress = progress.add_task("[green]Validating barcodes", total=linenum)
        # check for duplicates
        sort_out = subprocess.Popen(["sort", infile], stdout=subprocess.PIPE)
        dup_out = subprocess.run(["uniq", "-d"], stdin=sort_out.stdout, capture_output=True, text=True)
        if dup_out.stdout:
            print_error("duplicate barcodes", f"Duplicate barcodes were detected in {infile}, which will result in misleading simulated data.")
            print_solution_with_culprits("Check that you remove duplicate barcodes from your input file.", "Duplicates identified:")
            click.echo(dup_out.stdout, file = sys.stderr)
            sys.exit(1)
        for line,bc in enumerate(bc_file):
            length = validate(line + 1, bc)
            if length > limit:
                print_error("barcodes too long", f"Barcodes in [blue]{infile}[/blue] are [yellow]{length}bp[/yellow] and cannot exceed a length of [bold]{limit}bp[/bold]. Please use shorter barcodes.")
                sys.exit(1)
            lengths.add(length)
            progress.advance(task_progress)
    if len(lengths) > 1:
        print_error("incorrect format", f"Barcodes in [blue]{infile}[/blue] must all be a single length, but multiple lengths were detected: [yellow]" + ", ".join(lengths) + "[/yellow]")
        sys.exit(1)
    if return_len:
        return lengths.pop()
