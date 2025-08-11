"""Module with python-click types for command-line level validations of inputs"""

import os
import click
import pysam
import re
from harpy.common.misc import is_gzip
from harpy.common.printing import print_error, print_notice, print_solution, print_solution_offenders
import yaml
from pathlib import Path

class SAMfile(click.ParamType):
    """A CLI class to validate a BAM/SAM file as input. Checks for presence, format, and returns the absolute path"""
    name = "bam_file"
    def __init__(self, dir_ok: bool = True):
        super().__init__()
        self.dir_ok = dir_ok
    def convert(self, value, param, ctx):
        filepath = Path(value)
        if filepath.is_dir() and not self.dir_ok:
            self.fail("Alignment input cannot be a directory", param, ctx)
        if not filepath.is_dir():
            if not filepath.exists():
                self.fail(f"Alignment file {value} was not found.", param, ctx)
            if not os.access(value, os.R_OK):
                self.fail(f"Alignment file {value} does not have reading permission.", param, ctx)
                infiles = [filepath]
        else:
            re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)
            infiles = []
            for i in filepath.glob("*"):
                if i.is_file() and re_ext.search(i.name):
                    _file = i.resolve().as_posix()
                    if not os.access(_file, os.R_OK):
                        self.fail(f"Alignment file {_file} does not have reading permission.", param, ctx)
                    infiles.append(_file)
            if len(infiles) < 1:
                self.fail(f"There were no files ending with the accepted alignment extensions .bam/.sam (case insensitive) in {value}.")
        for i in infiles:
            try:
                with pysam.AlignmentFile(i, 'r', require_index=False):
                    pass
            except (ValueError, OSError):
                self.fail(f"File {value} was not recognized as being in BAM/SAM format.", param, ctx)
        return infiles

class FASTAfile(click.ParamType):
    """A CLI class to validate a FASTA file as input. Checks for presence, format, and returns the absolute path"""
    name = "fasta_file"
    def convert(self, value, param, ctx):
        filepath = Path(value)
        if not filepath.exists():
            self.fail(f"FASTA file {value} was not found.", param, ctx)
        if not filepath.is_file():
            self.fail(f"FASTA file {value} was not found.", param, ctx)
        if not os.access(value, os.R_OK):
            self.fail(f"FASTA file {value} does not have reading permission.", param, ctx)
        try:
            with pysam.FastxFile(value, persist=False) as f:
                for i in f:
                    if not i.name or i.quality:
                        raise ValueError
                    else:
                        return filepath.resolve().as_posix()
        except ValueError:
            self.fail(f"File {value} was not recognized as being in FASTA format.", param, ctx)
        
class FASTQfile(click.ParamType):
    """A CLI class to validate a FASTQ file as input. Checks for presence, format, and returns the absolute path"""
    name = "fastq_file"
    def __init__(self, dir_ok: bool = True):
        super().__init__()
        self.dir_ok = dir_ok

    def convert(self, value, param, ctx):
        filepath = Path(value)
        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        if filepath.is_dir() and not self.dir_ok:
            self.fail("FASTQ input cannot be a directory", param, ctx)
        if not filepath.is_dir():
            if not filepath.exists():
                self.fail(f"FASTQ file {value} was not found.", param, ctx)
            if not os.access(value, os.R_OK):
                self.fail(f"FASTQ file {value} does not have reading permission.", param, ctx)
            infiles = [filepath]
        else:
            infiles = []
            for i in filepath.glob("*"):
                if i.is_file() and re_ext.search(i.name):
                    _file = i.resolve().as_posix()
                    if not os.access(_file, os.R_OK):
                        self.fail(f"FASTQ file {_file} does not have reading permission.", param, ctx)
                    infiles.append(_file)
            if len(infiles) < 1:
                self.fail(f"There were no files ending with the accepted FASTQ extensions .fq[.gz] or .fastq[.gz] (case insensitive) in {value}.")

        for fastq in infiles:
            try:
                with pysam.FastxFile(fastq, persist=False) as f:
                    for i in f:
                        if not i.name or not i.quality:
                            raise ValueError
                        break
            except (ValueError, OSError):
                self.fail(f"File {fastq} was not recognized as being in FASTQ format.", param, ctx)
        return infiles

class VCFfile(click.ParamType):
    """A CLI class to validate a VCF/BCF file as input. Checks for presence, format, and returns the absolute path"""
    name = "vcf_file"
    def __init__(self, gzip_ok: bool = True):
        super().__init__()
        self.gzip_ok = gzip_ok

    def convert(self, value, param, ctx):
        filepath = Path(value)
        if not filepath.exists():
            self.fail(f"Variant call format file {value} was not found", param, ctx)
        if not os.access(value, os.R_OK):
            self.fail(f"Variant call format file {value} does not have read permission.", param, ctx)
        if is_gzip(value) and not self.gzip_ok:
            self.fail(f"Variant call format file {value} cannot be gzipped. Please decompress it.", param, ctx)
        try:
            with pysam.VariantFile(value, 'r'):
                return filepath.resolve().as_posix()
        except (ValueError, OSError):
            self.fail(f"File {value} was not recognized as being in BCF/VCF[.GZ] format.", param, ctx)

class PopulationFile(click.ParamType):
    name = "populations_file"

    def convert(self, value, param, ctx):
        filepath = Path(value)
        if not filepath.exists():
            self.fail(f"Sample grouping file {value} was not found", param, ctx)
        if not os.access(value, os.R_OK):
            self.fail(f"Sample grouping file {value} does not have read permission.", param, ctx)

        with open(value, "r", encoding="utf-8") as f:
            rows = [i for i in f.readlines() if i != "\n" and not i.lstrip().startswith("#")]
            invalids = [(i,j) for i,j in enumerate(rows) if len(j.split()) < 2]
            if invalids:
                print_error(
                    "invalid format",
                    f"There are [bold]{len(invalids)}[/] rows in [bold]{value}[/] without a space/tab delimiter or don't have two entries for sample[dim]<tab>[/dim]population. Terminating Harpy to avoid downstream errors.",
                    False
                )
                print_solution_offenders(
                    f"Make sure every entry in [bold]{value}[/] uses space or tab delimeters and has both a sample name and population designation. You may comment out rows with a [green bold]#[/] to have Harpy ignore them.",
                    "The rows and values causing this error are",
                    [f"{i[0]+1}\t{i[1]}" for i in invalids]
                )
        return filepath.resolve().as_posix()

class InputFile(click.ParamType):
    """A class for a click type that verifies that a file exists and that it has an expected extension. Returns the absolute path"""
    name = "input_file"
    def __init__(self, filetype, gzip_ok):
        super().__init__()
        self.filetype = filetype
        self.gzip_ok = gzip_ok
    def convert(self, value, param, ctx):
        filedict = {
            "fasta": [".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".frn"],
            "vcf": ["vcf", "bcf", "vcf.gz"],
            "gff": [".gff",".gff3"]
        }
        if self.filetype not in filedict:
            self.fail(f"Extension validation for {self.filetype} is not yet implemented. This error should only appear during development; if you are a user and seeing this, please post an issue on GitHub: https://github.com/pdimens/harpy/issues/new?assignees=&labels=bug&projects=&template=bug_report.yml")
        if not os.path.exists(value):
            self.fail(f"{value} does not exist. Please check the spelling and try again.", param, ctx)
        elif not os.access(value, os.R_OK):
            self.fail(f"{value} is not readable. Please check file/directory permissions and try again", param, ctx)
        if os.path.isdir(value):
            self.fail(f"{value} is a directory, but input should be a file.", param, ctx)
        valid = False
        lowercase = value.lower()
        for ext in filedict[self.filetype]:
            valid = True if lowercase.endswith(ext) else valid
            if self.gzip_ok:
                valid = True if lowercase.endswith(ext + ".gz") else valid
        if not valid and not self.gzip_ok:
                self.fail(f"{value} does not end with one of the expected extensions [" + ", ".join(filedict[self.filetype]) + "]. Please verify this is the correct file type and rename the extension for compatibility.", param, ctx)
        if not valid and self.gzip_ok:
            self.fail(f"{value} does not end with one of the expected extensions [" + ", ".join(filedict[self.filetype]) + "]. Please verify this is the correct file type and rename the extension for compatibility. Gzip compression (ending in .gz) is allowed.", param, ctx)
        return Path(value).resolve().as_posix()

class HPCProfile(click.ParamType):
    """A class for a click type which accepts a file with a snakemake HPC profile. Does validations to make sure it's the config file and not the directory."""
    name = "hpc_profile"
    def convert(self, value, param, ctx):
        if not os.path.exists(value):
            self.fail(f"{value} does not exist. Please check the spelling and try again.", param, ctx)
        if os.path.isdir(value):
            self.fail(f"{value} is a directory, but input should be a yaml file.", param, ctx)
        if not os.access(value, os.R_OK):
            self.fail(f"{value} is not readable. Please check file permissions and try again", param, ctx)
        with open(value, "r") as file:
            try:
                yaml.safe_load(file)
            except yaml.YAMLError as exc:
                self.fail(f"Formatting error in {value}: {exc}")
        return Path(value).resolve().as_posix()

class DemuxSchema(click.ParamType):
    """A class for a click type that accepts a demultiplex schema and performs validation"""
    name = "demultiplex_schema"

    def convert(self, value, param, ctx):
        if not os.path.exists(value):
            self.fail(f"{value} does not exist. Please check the spelling and try again.", param, ctx)
        if os.path.isdir(value):
            self.fail(f"{value} is a directory, but should be a file.", param, ctx)
        if not os.access(value, os.R_OK):
            self.fail(f"{value} is not readable. Please check file permissions and try again", param, ctx)

        # can be Axx, Bxx, Cxx, Dxx
        code_letters = set()
        segment_ids = set()
        samples = set()
        duplicates = False
        segment_pattern = re.compile(r'^[A-D]\d{2}$')
        with open(value, 'r') as file:
            for line in file:
                try:
                    sample, segment_id = line.rstrip().split()
                    if not segment_pattern.match(segment_id):
                        print_error("invalid segment format", f"Segment ID [green]{segment_id}[/] does not follow the expected format.", False)
                        print_solution("This haplotagging design expects segments to follow the format of letter [green bold]A-D[/] followed by [bold]two digits[/], e.g. [green bold]C51[/]). Check that your ID segments or formatted correctly and that you are attempting to demultiplex with a workflow appropriate for your data design.")
                    code_letters.add(segment_id[0])
                    if sample in samples:
                        duplicates = True
                    samples.add(sample)
                    if segment_id in segment_ids:
                        print_error("ambiguous segment ID", "An ID segment must only be associated with a single sample.", False)
                        print_solution_offenders(
                            "A barcode segment can only be associated with a single sample. For example: [green bold]C05[/] cannot identify both [green]sample_01[/] and [green]sample_2[/]. In other words, a segment can only appear once.",
                            "The segment triggering this error is",
                            segment_id
                        )
                    else:
                        segment_ids.add(segment_id)
                except ValueError:
                    # skip rows without two columns
                    continue
        if not code_letters:
            print_error("incorrect schema format", f"Schema file [blue]{os.path.basename(value)}[/] has no valid rows. Rows should be sample<tab>segment, e.g. sample_01<tab>C75")
        if len(code_letters) > 1:
            print_error("invalid schema", f"Schema file [blue]{os.path.basename(value)}[/] has sample IDs occurring in different barcode segments.", False)
            print_solution_offenders(
                "All sample IDs for this barcode design should be in a single segment, such as [bold green]C[/] or [bold green]D[/]. Make sure the schema contains only one segment.",
                "The segments identified in the schema",
                ", ".join(code_letters)
            )
        if duplicates:
            print_notice("Sample names appear more than once, assuming this was intentional")
        return Path(value).resolve().as_posix()
