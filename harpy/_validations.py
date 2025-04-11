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

# not currently used, but keeping it here for posterity
def is_bgzipped(file_path: str) -> bool:
    """Check if a file is truly BGZF-compressed (not just GZIP) by looking for the BGZF EOF marker."""
    try:
        # Try reading the BGZF EOF marker (last 28 bytes)
        with open(file_path, 'rb') as f:
            f.seek(-28, 2)  # Seek to 28 bytes before end
            eof_block = f.read()
            # BGZF EOF marker signature
            bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
            return eof_block == bgzf_eof
    except (IOError, OSError):
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
                culprit_text.append(f"[red]Invalid columns:[/] " + ", ".join(invalid_col))
            if missing_col:
                errtext.append("\n  - missing columns")
                culprit_text.append( f"[yellow]Missing columns:[/] " + ", ".join(missing_col))

            print_error("incorrect columns", f"Parameter file [bold]{parameters}[/] has incorrect column names" + "".join(errtext) + f"\nValid names are: [green bold]" + " ".join(colnames) + "[/]")
            print_solution_with_culprits(
                f"Fix the headers in [bold]{parameters}[/] or use [blue bold]harpy template[/] to generate a valid parameter file and modify it with appropriate values.",
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
            print_error("invalid rows", f"Parameter file [blue]{parameters}[/] is formatted incorrectly. Not all rows have the expected {n_cols} columns.")
            print_solution_with_culprits(
                f"See the problematic rows below. Check that you are using a whitespace (space or tab) delimeter in [blue]{parameters}[/] or use [blue green]harpy template[/blue green] to generate a valid parameter file and modify it with appropriate values.",
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
            print_error("invalid parameter values", f"Parameter file [bold]{parameters}[/] is formatted incorrectly. Some rows have incorrect values for one or more parameters.")
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
        print_error("sample ID mismatch",f"There are [bold]{len(culpritfiles)}[/] alignment files whose RG tags do not match their filenames.")
        print_solution_with_culprits(
            "For alignment files (sam/bam), the base of the file name must be identical to the [green bold]@RD ID:[/] tag in the file header. For example, a file named \'sample_001.bam\' should have the [green bold]@RG ID:sample_001[/] tag. Use the [blue bold]renamebam[/] script to properly rename alignment files so as to also update the @RG header.",
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
        print_error("vcf not phased", "The input variant file needs to be phased into haplotypes, but no [green]FORMAT/PS[/] or [green]FORMAT/HP[/] fields were found.")
        print_solution(f"Phase [bold]{bn}[/] into haplotypes using [blue bold]harpy phase[/] or another manner of your choosing and use the phased vcf file as input. If you are confident this file is phased, then the phasing does not follow standard convention and you will need to make sure the phasing information appears as either [green]FORMAT/PS[/] or [green]FORMAT/HP[/] tags.")
        sys.exit(1)

def validate_popfile(infile: str) -> None:
    """Validate the input population file to make sure there are two entries per row"""
    with open(infile, "r", encoding="utf-8") as f:
        rows = [i for i in f.readlines() if i != "\n" and not i.lstrip().startswith("#")]
        invalids = [(i,j) for i,j in enumerate(rows) if len(j.split()) < 2]
        if invalids:
            print_error("invalid format", f"There are [bold]{len(invalids)}[/] rows in [bold]{infile}[/] without a space/tab delimiter or don't have two entries for sample[dim]<tab>[/dim]population. Terminating Harpy to avoid downstream errors.")
            print_solution_with_culprits(
                f"Make sure every entry in [bold]{infile}[/] uses space or tab delimeters and has both a sample name and population designation. You may comment out rows with a [green bold]#[/] to have Harpy ignore them.",
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
        print_error("mismatched inputs", f"There are [bold]{len(missing_samples)}[/] samples found in [blue]{fromthis}[/] that are not in [blue]{inthis}[/]. Terminating Harpy to avoid downstream errors.")
        print_solution_with_culprits(
            f"[blue]{fromthis}[/] cannot contain samples that are absent in [blue]{inthis}[/]. Check the spelling or remove those samples from [blue]{fromthis}[/] or remake the vcf file to include/omit these samples. Alternatively, toggle [green]--vcf-samples[/] to aggregate the sample list from the input files or [blue]{vcf}[/].",
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
        print_error("contigs absent", f"Some of the provided contigs were not found in [blue]{shortname}[/]. This will definitely cause plotting errors in the workflow.")
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
        print_notice(f"There are [bold]{len(overlooked)}[/] samples found in the inputs that weren\'t included in [blue]{popfile}[/]. This will [bold]not[/] cause errors and can be ignored if it was deliberate. Commenting or removing these lines will avoid this message. The samples are:\n" + ", ".join(overlooked))
    if len(missing_samples) > 0:
        print_error("mismatched inputs", f"There are [bold]{len(missing_samples)}[/] samples included in [blue]{popfile}[/] that weren\'t found in in the inputs. Terminating Harpy to avoid downstream errors.")
        print_solution_with_culprits(
            f"Make sure the spelling of these samples is identical in the inputs and [blue]{popfile}[/], or remove them from [blue]{popfile}[/].",
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
                    print_error("invalid segment format", f"Segment ID [green]{segment_id}[/] does not follow the expected format.")
                    print_solution("This haplotagging design expects segments to follow the format of letter [green bold]A-D[/] followed by [bold]two digits[/], e.g. [green bold]C51[/]). Check that your ID segments or formatted correctly and that you are attempting to demultiplex with a workflow appropriate for your data design.")
                    sys.exit(1)
                code_letters.add(segment_id[0])
                samples.add(sample)
                if segment_id in segment_ids:
                    print_error("ambiguous segment ID", "An ID segment must only be associated with a single sample.")
                    print_solution_with_culprits(
                        "A barcode segment can only be associated with a single sample. For example: [green bold]C05[/] cannot identify both [green]sample_01[/] and [green]sample_2[/]. In other words, a segment can only appear once.",
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
        print_error("incorrect schema format", f"Schema file [blue]{os.path.basename(infile)}[/] has no valid rows. Rows should be sample<tab>segment, e.g. sample_01<tab>C75")
        sys.exit(1)
    if len(code_letters) > 1:
        print_error("invalid schema", f"Schema file [blue]{os.path.basename(infile)}[/] has sample IDs occurring  in different barcode segments.")
        print_solution_with_culprits(
            "All sample IDs for this barcode design should be in a single segment, such as [bold green]C[/] or [bold green]D[/]. Make sure the schema contains only one segment.",
            "The segments identified in the schema:"
        )
        click.echo(", ".join(code_letters))
        sys.exit(1)
    if return_len:
        return len(samples)

def validate_regions(regioninput: int | str, genome: str) -> None:
    """validates the --regions input of harpy snp to infer whether it's an integer, region, or file"""
    def contig_lens(gen):
        # read in contigs for length/presence checks
        contigs = {}
        with safe_read(gen) as fopen:
            # get contig lengths
            for line in fopen:
                if line.startswith(">"):
                    cn = line.rstrip("\n").lstrip(">").split()[0]
                    contigs[cn] = 0
                else:
                    contigs[cn] += len(line.rstrip("\n"))
        return contigs

    try:
        # is an int
        int(regioninput)
        return
    except ValueError:
        pass

    # is a file specifying regions
    if os.path.isfile(regioninput):
        contigs = contig_lens(genome)
        with open(regioninput, "r", encoding="utf-8") as fin:
            for idx, line in enumerate(fin, 1):
                row = line.split()
                if len(row) != 3:
                    print_error("invalid format", f"The input file is formatted incorrectly at line {idx}. This is the first row triggering this error, but it may not be the only one.")
                    print_solution_with_culprits(
                        f"Rows in [blue]{regioninput}[/] need to be [bold]space[/] or [bold]tab[/] delimited with the format [yellow bold]contig start end[/] where [yellow bold]start[/] and [yellow bold]end[/] are integers.",
                        "Rows triggering this error:"
                        )
                    click.echo(line, file = sys.stderr)
                    sys.exit(1)
                else:
                    try:
                        start = int(row[1])
                        end = int(row[2])
                    except ValueError:
                        print_solution_with_culprits(
                            f"Rows in [blue]{regioninput}[/] need to be [bold]space[/] or [bold]tab[/] delimited with the format [yellow bold]contig start end[/] where [yellow bold]start[/] and [yellow bold]end[/] are integers.",
                            "Rows triggering this error:"
                            )
                        click.echo(line, file = sys.stderr)
                        sys.exit(1)
                if row[0] not in contigs:
                    print_error("missing contig", f"The contig listed at row {idx} ([bold yellow]{row[0]}[/]) is not present in ([blue]{os.path.basename(genome)}[/]). This is the first row triggering this error, but it may not be the only one.")
                    print_solution(
                        f"Check that all the contigs listed in [blue]{os.path.basename(regioninput)}[/] are also present in [blue]{os.path.basename(genome)}[/]",
                        "Row triggering this error:"
                    )
                    click.echo(line, file = sys.stderr)
                    sys.exit(1)
                if start > end:
                    print_error("invalid interval", f"The interval start position is greater than the interval end position at row {idx}. This is the first row triggering this error, but it may not be the only one.")
                    print_solution(
                        f"Check that all rows in [blue]{os.path.basename(regioninput)}[/] have a [bold yellow]start[/] position that is less than the [bold yellow]end[/] position.",
                        "Row triggering this error:"
                    )
                    click.echo(line, file = sys.stderr)
                    sys.exit(1)
                if start > contigs[row[0]] or end > contigs[row[0]]:
                    print_error("invalid interval", f"The interval start or end position is out of bounds at row {idx}. This is the first row triggering this error, but it may not be the only one.")
                    print_solution(
                        f"Check that the intervals present in [blue]{os.path.basename(regioninput)}[/] are within the bounds of the lengths of their respective contigs. This specific error is triggered for [bold yellow]{row[0]}[/], which has a total length of [bold]{contigs[row[0]]}[/].",
                        "Row triggering this error:"
                    )
                    click.echo(line, file = sys.stderr)
                    sys.exit(1)
        return

    contig,positions = regioninput.split(":")
    startpos,endpos = [int(i) for i in positions.split("-")]
    contigs = contig_lens(genome)
    if contig not in contigs:
        print_error("contig not found", f"The contig ([bold yellow]{contig}[/]) of the input region [yellow bold]{regioninput}[/] was not found in [blue]{genome}[/].")
        sys.exit(1)
    if startpos > contigs[contig]:
        print_error("region out of bounds", f"The start position ([yellow bold]{startpos}[/]) exceeds the total length of contig [yellow bold]{contig}[/] ({contigs[contig]})")
        sys.exit(1)
    if endpos > contigs[contig]:
        print_error("region out of bounds", f"The end position ([yellow bold]{endpos}[/]) exceeds the total length of contig [yellow bold]{contig}[/] ({contigs[contig]})")
        sys.exit(1)
    return

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
                    print_error("consecutive contig names", f"All contig names must be followed by at least one line of nucleotide sequences, but two consecutive lines of contig names were detected. This issue was identified at line [bold]{line_num}[/] in [blue]{genofile}[/], but there may be others further in the file.")
                    print_solution("See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
                    sys.exit(1)
                else:
                    last_header = True
                if len(line.rstrip()) == 1:
                    print_error("unnamed contigs", f"All contigs must have an alphanumeric name, but a contig was detected without a name. This issue was identified at line [bold]{line_num}[/] in [blue]{genofile}[/], but there may be others further in the file.")
                    print_solution("See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
                    sys.exit(1)
                if line.startswith("> "):
                    print_error("invalid contig names", f"All contig names must be named [green bold]>contig_name[/], without a space, but a contig was detected with a space between the [green bold]>[/] and contig_name. This issue was identified at line [bold]{line_num}[/] in [blue]{genofile}[/], but there may be others further in the file.")
                    print_solution("See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
                    sys.exit(1)
            elif line == "\n":
                print_error("empty lines", f"Empty lines are not permitted in FASTA files, but one was detected at line [bold]{line_num}[/] in [blue]{genofile}[/]. The scan ended at this error, but there may be others further in the file.")
                print_solution("See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
                sys.exit(1)
            else:
                seq += 1
                last_header = False
    solutiontext = "FASTA files must have at least one contig name followed by sequence data on the next line. Example:\n"
    solutiontext += "[green]  >contig_name\n  ATACAGGAGATTAGGCA[/]\n"
    # make sure there is at least one of each
    if seq_id == 0:
        print_error("contig names absent", f"No contig names detected in [blue]{genofile}[/].")
        print_solution(solutiontext + "\nSee the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/")
        sys.exit(1)
    if seq == 0:
        print_error("sequences absent", f"No sequences detected in [blue]{genofile}[/].")
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
        print_error("contigs absent", f"Some of the provided contigs were not found in [blue]{shortname}[/]. This will definitely cause plotting errors in the workflow.")
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
                    print_error("clashing barcode tags", f"Both [green bold]BC:Z[/] and [green bold]BX:Z[/] tags were detected in the read headers for [blue]{os.path.basename(fastq)}[/]. Athena accepts [bold]only[/] one of [green bold]BC:Z[/] or [green bold]BX:Z[/].")
                    print_solution("Check why your data has both tags in use and remove/rename one of the tags.")
                    sys.exit(1)
            # check for one or the other after parsing is done
            if not BX and not BC:
                print_error("no barcodes found",f"No [green bold]BC:Z[/] or [green bold]BX:Z[/] tags were detected in read headers for [blue]{os.path.basename(fastq)}[/]. Athena requires the linked-read barcode to be present as either [green bold]BC:Z[/] or [/]BX:Z[/] tags.")
                print_solution("Check that this is linked-read data and that the barcode is demultiplexed from the sequence line into the read header as either a `BX:Z` or `BC:Z` tag.")
                sys.exit(1)

    # parellelize over the fastq list
    with harpy_progressbar(quiet) as progress, ThreadPoolExecutor(max_workers=threads) as executor:
        task_progress = progress.add_task("[green]Validating FASTQ inputs", total=len(fastq_list))
        futures = [executor.submit(validate, i) for i in fastq_list]
        for future in as_completed(futures):
            progress.advance(task_progress)

def validate_barcodefile(infile: str, return_len: bool = False, quiet: int = 0, limit: int = 60, gzip_ok: bool = True, haplotag_only: bool = False, check_dups: bool = True) -> None | int:
    """Does validations to make sure it's one length, within a length limit, one per line, and nucleotides"""
    if is_gzip(infile) and not gzip_ok:
        print_error("incorrect format", f"The input file must be in uncompressed format. Please decompress [blue]{infile}[/] and try again.")
        sys.exit(1)
    lengths = set()
    nucleotides = {'A','C','G','T'}
    def validate(line_num, bc_line):
        barcode = bc_line.rstrip()
        if len(barcode.split()) > 1:
            print_error("incorrect format", f"There must be one barcode per line, but multiple entries were detected on [bold]line {line_num}[/] in [blue]{infile}[/]")
            sys.exit(1)
        if not set(barcode).issubset(nucleotides) or barcode != barcode.upper():
            print_error("incorrect format", f"Invalid barcode format on [bold]line {line_num }[/]: [yellow]{barcode}[/].\nBarcodes in [blue]{infile}[/] must be captial letters and only contain standard nucleotide characters [green]ATCG[/].")
            sys.exit(1)
        return len(barcode)

    with safe_read(infile) as bc_file, harpy_progressbar(quiet) as progress:
        out = subprocess.Popen(['wc', '-l', infile], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
        linenum = int(out.partition(b' ')[0])
        if linenum > 96**4 and haplotag_only:
            print_error("Too many barcodes", f"The maximum number of barcodes possible with haplotagging is [bold]96^4 (84,934,656)[/], but there are [yellow]{linenum}[/] barcodes in [blue]{infile}[/]. Please use fewer barcodes.")
            sys.exit(1)
        task_progress = progress.add_task("[green]Validating barcodes", total=linenum)
        # check for duplicates
        if check_dups:
            sort_out = subprocess.Popen(["sort", infile], stdout=subprocess.PIPE)
            dup_out = subprocess.run(["uniq", "-d"], stdin=sort_out.stdout, capture_output=True, text=True)
            if dup_out.stdout:
                print_error("duplicate barcodes", f"Duplicate barcodes were detected in {infile}, which will result in misleading simulated data.")
                print_solution_with_culprits("Check that you remove duplicate barcodes from your input file.", "Duplicates identified:")
                click.echo(dup_out.stdout, file = sys.stderr)
                sys.exit(1)
        for line,bc in enumerate(bc_file, 1):
            length = validate(line, bc)
            if length > limit:
                print_error("barcodes too long", f"Barcodes in [blue]{infile}[/] are [yellow]{length}bp[/] and cannot exceed a length of [bold]{limit}bp[/]. Please use shorter barcodes.")
                sys.exit(1)
            lengths.add(length)
            if len(lengths) > 1:
                print_error("inconsistent length", f"Barcodes in [blue]{infile}[/] must all be a single length, but multiple lengths were detected: [yellow]" + ", ".join(lengths) + "[/]")
                sys.exit(1)
            progress.advance(task_progress)
    if return_len:
        return lengths.pop()

