"""Module with python-click types for command-line level validations of inputs"""

import click
import os
from pathlib import Path
import re
import yaml
from harpy.common.file_ops import is_gzip
from harpy.common.printing import HarpyPrint

hp = HarpyPrint()

class SAMfile(click.ParamType):
    """A CLI class to validate a BAM/SAM file as input. Checks for presence, format, and returns the absolute path"""
    name = "bam_file"
    def __init__(self, dir_ok: bool = True, single: bool = False):
        super().__init__()
        self.single = single
        self.dir_ok = False if single else dir_ok
        self.re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)
        self.inv_pattern = re.compile(r'[^a-z0-9._-]+', re.IGNORECASE)

    def convert(self, value, param, ctx):
        infiles = []
        filepath = Path(value)
        if not filepath.exists():
            self.fail(f"{value} was not found.", param, ctx)

        if not filepath.is_dir():
            _f = filepath.resolve().as_posix()
            if not self.re_ext.search(_f):
                self.fail(f"{value} does not end with the accepted extensions for alignment files: .bam/.sam (case insensitive).")
            infiles.append(_f)
        else:
            if filepath.is_dir() and not self.dir_ok:
                self.fail("Alignment input cannot be a directory", param, ctx)
            for i in filepath.glob("*"):
                if i.is_file() and self.re_ext.search(i.name):
                    infiles.append(i.resolve().as_posix())

        # name and permission validations
        for _file in infiles:
            if not os.access(_file, os.R_OK):
                self.fail(f"Alignment file {_file} does not have user read permission.", param, ctx)
            #if self.inv_pattern.search(_file):
            if self.inv_pattern.search(os.path.basename(_file)):
                self.fail(f"Invalid characters detected in file name {_file}. Valid file names may contain only:\n  - A-Z 0-9 characters (case insensitive)\n  - . (period)\n  - _ (underscore)\n  - - (dash)")

        if len(infiles) < 1:
            self.fail(f"There were no files ending with the accepted alignment extensions .bam/.sam (case insensitive) in {value}.")

        if self.single:
            return infiles[0]
        else:
            return infiles

class FASTAfile(click.ParamType):
    """A CLI class to validate a FASTA file as input. Checks for presence, format, and returns the absolute path"""
    name = "fasta_file"

    def convert(self, value, param, ctx):
        re_ext = re.compile(r"\.(fa|fas|fasta|fna|ffn|frn)(?:\.gz)?$", re.IGNORECASE)
        filepath = Path(value)
        if not filepath.exists() or not filepath.is_file():
            self.fail(f"FASTA file {value} was not found.", param, ctx)
        if not os.access(value, os.R_OK):
            self.fail(f"FASTA file {value} does not have reading permission.", param, ctx)

        if not re_ext.search(value):
            self.fail(f"File {value} does not have any of the recognized [case insensitive] FASTA file extensions (.fa, .fas, .fasta, .fna, .ffn, .frn). Gzipping (.gz) is also permitted.", param, ctx)

        return filepath.resolve().as_posix()
        
class FASTQfile(click.ParamType):
    """
    A CLI class to validate a FASTQ file as input. Checks for presence, format, and returns the absolute path.
    Setting single to True returns a str, otherwise returns a list.
    """
    name = "fastq_file"
    def __init__(self, dir_ok: bool = True, single: bool = False):
        super().__init__()
        self.single = single
        self.dir_ok = False if single else dir_ok
        self.inv_pattern = re.compile(r'[^a-z0-9._-]+', re.IGNORECASE)
        self.re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)

    def convert(self, value, param, ctx):
        infiles = []
        filepath = Path(value)
        if not filepath.exists():
            self.fail(f"{value} was not found.", param, ctx)
        if filepath.is_dir() and not self.dir_ok:
            self.fail("FASTQ input cannot be a directory", param, ctx)

        if not filepath.is_dir():
            _f = filepath.resolve().as_posix()
            if not self.re_ext.search(_f):
                self.fail(f"{value} does not end with the accepted extensions for FASTQ files: .fq[.gz]/.fastq[.gz] (case insensitive).")
            infiles.append(filepath.resolve().as_posix())
        else:
            for i in filepath.glob("*"):
                if i.is_file() and self.re_ext.search(i.name):
                    infiles.append(i.resolve().as_posix())

        for _file in infiles:
            if not os.access(_file, os.R_OK):
                self.fail(f"FASTQ file {_file} does not have user read permission.", param, ctx)
            if self.inv_pattern.search(os.path.basename(_file)):
                self.fail(f"Invalid characters detected in file name {_file}. Valid file names may contain only:\n  - A-Z 0-9 characters (case insensitive)\n  - . (period)\n  - _ (underscore)\n  - - (dash)")

        if len(infiles) < 1:
            self.fail(f"There were no files ending with the accepted FASTQ extensions .fq[.gz] or .fastq[.gz] (case insensitive) in {value}.")

        if self.single:
            return infiles[0]
        else:
            return infiles

class VCFfile(click.ParamType):
    """A CLI class to validate a VCF/BCF file as input. Checks for presence, format, and returns the absolute path"""
    name = "vcf_file"
    def __init__(self, gzip_ok: bool = True):
        super().__init__()
        self.gzip_ok = gzip_ok
        self.re_ext = re.compile(r"\.(vcf|vcf\.gz|bcf)$", re.IGNORECASE)

    def convert(self, value, param, ctx):
        filepath = Path(value)
        _file = filepath.resolve().as_posix()
        if not filepath.exists():
            self.fail(f"Variant call format file {value} was not found", param, ctx)
        if not os.access(value, os.R_OK):
            self.fail(f"Variant call format file {value} does not have read permission.", param, ctx)
        if is_gzip(value) and not self.gzip_ok:
            self.fail(f"Variant call format file {value} cannot be gzipped. Please decompress it.", param, ctx)

        if not self.re_ext.search(_file):
            self.fail(f"File {value} was not recognized as having the [case-insensitive] accepted variant call format file extensions: .vcf, .vcf.gz, .bcf", param, ctx)
        else:
            return _file

class PopulationFile(click.ParamType):
    name = "populations_file"

    def convert(self, value, param, ctx):
        filepath = Path(value)
        if not filepath.exists():
            self.fail(f"Sample grouping file {value} was not found", param, ctx)
        if not os.access(value, os.R_OK):
            self.fail(f"Sample grouping file {value} does not have read permission.", param, ctx)
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
        filepath = Path(value)
        if not filepath.exists():
            self.fail(f"{value} does not exist. Please check the spelling and try again.", param, ctx)
        if filepath.is_dir():
            self.fail(f"{value} is a directory, but should be a file.", param, ctx)
        if not os.access(value, os.R_OK):
            self.fail(f"{value} is not readable. Please check file permissions and try again", param, ctx)

        # can be Axx, Bxx, Cxx, Dxx
        code_letters = set()
        segment_ids = set()
        samples = set()
        duplicates = False
        segment_pattern = re.compile(r'^[A-D]\d{2}$')
        with open(filepath, 'r') as file:
            for line in file:
                try:
                    sample, segment_id = line.rstrip().split()
                    if not segment_pattern.match(segment_id):
                        hp.error(
                            "invalid segment format",
                            f"Segment ID [green]{segment_id}[/] does not follow the expected format.",
                            "This haplotagging design expects segments to follow the format of letter [green bold]A-D[/] followed by [bold]two digits[/], e.g. [green bold]C51[/]). Check that your ID segments or formatted correctly and that you are attempting to demultiplex with a workflow appropriate for your data design."
                        )
                    code_letters.add(segment_id[0])
                    if sample in samples:
                        duplicates = True
                    samples.add(sample)
                    if segment_id in segment_ids:
                        hp.error(
                            "ambiguous segment ID",
                            "An ID segment must only be associated with a single sample.",                        
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
            hp.error("incorrect schema format", f"Schema file [blue]{os.path.basename(value)}[/] has no valid rows. Rows should be sample<tab>segment, e.g. sample_01<tab>C75")
        if len(code_letters) > 1:
            hp.error(
                "invalid schema",
                f"Schema file [blue]{os.path.basename(value)}[/] has sample IDs occurring in different barcode segments.",
                "All sample IDs for this barcode design should be in a single segment, such as [bold green]C[/] or [bold green]D[/]. Make sure the schema contains only one segment.",
                "The segments identified in the schema",
                ", ".join(code_letters)
            )
        if duplicates:
            hp.notice("Sample names appear more than once, assuming this was intentional")
        return filepath.resolve().as_posix()

class impute_strategy(click.ParamType):
    name = "impute_strategy"

    def isInt(self, x: str) -> bool:
        try:
            int(x)
        except ValueError:
            return False
        return True

    def convert(self, value, param, ctx):
        if value == "all":
            return value

        pieces = value.split(":")
        if len(pieces) != 2:
            self.fail(f"{value} is not in a recognized format. Expected one of:\n  - all\n  - chrom:start-end (e.g., chr1:1-50000)\n  - window:size (e.g., window:1000000)", param, ctx)

        if pieces[0] == "window":
            try:
                _win = int(pieces[1])
                if _win < 100:
                    self.fail(f"Window size must be >= 100, but got {value}", param, ctx)  
            except ValueError:
                self.fail(f"The window strategy format is incorrect.\n  - expected window:size (e.g., window:1000000)\n  - got {value}")
        else:
            parts = pieces[1].split("-")
            if len(parts) != 2 or any([not self.isInt(i) for i in parts]):
                self.fail(f"Region must be contig:start-end (e.g., chr1:300-80000), but got {value}")
            start, end = map(int, parts)  
            if start < 1 or end < start:  
                self.fail(f"Region coordinates must satisfy start >= 1 and end >= start, but got {value}", param, ctx) 

        return value
