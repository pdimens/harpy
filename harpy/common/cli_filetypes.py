"""Module with python-click types for command-line level validations of inputs"""

import os
import click
import pysam
import re
from harpy.common.validations import is_gzip
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
                    infiles.append(i)
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
                    infiles.append(i)
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
        return value
