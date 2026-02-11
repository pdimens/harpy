"""Module with helper function to set up Harpy workflows"""

import glob
import gzip
import importlib.resources as resources
import os
from pathlib import Path
import pysam
import re
import shutil
import sys
from harpy.common.printing import HarpyPrint

def filepath(infile: str) -> str:
    """returns a posix-formatted absolute path of infile"""
    return Path(infile).resolve().as_posix()

def symlink(original: str, destination: str) -> None:
    """Create a symbolic link from original -> destination if the destination doesn't already exist."""
    if not (Path(destination).is_symlink() or Path(destination).exists()):
        Path(destination).symlink_to(Path(original).resolve())

def fetch_template(target: str, outfile = None) -> None:
    """
    Retrieve the target file from harpy.templates and print to outfile. Prints
    to stdout if no outfile provided.
    """
    source_file = resources.files("harpy.templates") / target
    if outfile:
        _dir = os.path.dirname(outfile)
        if _dir:
            os.makedirs(_dir, exist_ok=True)
        _out = open(outfile, "w")
    else:
        _out = sys.stdout
    try:
        if target.lower().endswith(".png"):
            shutil.copy(str(source_file), outfile)
        else:
            with resources.as_file(source_file) as _source, open(_source, 'r') as f:
                _out.write(f.read() + "\n")
    except (FileNotFoundError, KeyError):
        HarpyPrint().error(
            "template file missing",
            f"The required template file [blue bold]{target}[/] was not found within the Harpy installation.",
            "There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration."
        )
    finally:
        if outfile:
            _out.close()

def gzip_file(infile: str) -> None:
    """gzip a file and delete the original, using only python"""
    if os.path.exists(infile):
        with open(infile, 'rb') as f_in, gzip.open(infile + '.gz', 'wb', 6) as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(infile)

def last_sm_log(directory) -> str:
    '''Find and return the name of the last snakemake log file, otherwise return an empty string'''
    err_dir = os.path.join(directory, ".snakemake", "log")
    if not os.path.exists(err_dir) or not os.path.exists(directory):
        return ""
    files = [i for i in glob.iglob(f"{err_dir}/*.log*")]        
    return os.path.basename(sorted(files, key = os.path.getmtime)[-1])

def purge_empty_logs(output_directory):
    """scan target_dir and remove empty files, then scan it again and remove empty directories"""
    for logfile in glob.glob(f"{output_directory}/logs/**/*", recursive = True):
        if os.path.isfile(logfile) and os.path.getsize(logfile) == 0:
            os.remove(logfile)
    for logfile in glob.glob(f"{output_directory}/logs/**/*", recursive = True):
        if os.path.isdir(logfile) and not os.listdir(logfile):
            os.rmdir(logfile)

def safe_read(file_path: str):
    """returns the proper file opener for reading if a file_path is gzipped"""
    try:
        with gzip.open(file_path, 'rt') as f:
            f.read(10)
        return gzip.open(file_path, 'rt')
    except gzip.BadGzipFile:
        return open(file_path, 'r')

def is_gzip(file_path: str) -> bool:
    """helper function to determine if a file is gzipped"""
    try:
        with gzip.open(file_path, 'rt') as f:
            f.read(10)
        return True
    except (gzip.BadGzipFile, UnicodeDecodeError):
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

def pop_manifest(groupingfile, filelist) -> dict:
    '''
    Create dictionary of population => filenames to make it easier to
    set the snakemake rules/wildcards for population grouping
    '''
    d = {}
    with open(groupingfile) as f:
        for line in f:
            samp, pop = line.rstrip().split()
            if samp.lstrip().startswith("#"):
                continue
            r = re.compile(fr".*/({samp.lstrip()})\.(bam|sam)$", flags = re.IGNORECASE)
            sampl = list(filter(r.match, filelist))[0]
            if pop not in d.keys():
                d[pop] = [sampl]
            else:
                d[pop].append(sampl)
    return d

def naibr_extra(argsDict: dict, extra) -> dict:
    if extra:
        words = [i for i in re.split(r"\s|=", extra) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            if "blacklist" in i or "candidates" in i:
                argsDict[i[0].lstrip("-")] = i[1]
    return argsDict

def genomic_windows(input: str, output: str, window: int = 10000, mode: int = 1):
    """
    Create a BED file of fixed intervals from a fasta or fai file (generated with samtools faidx).
    Nearly identical to bedtools makewindows, except the intervals are nonoverlapping. `Window` is
    the interval size, `mode` is whether to make `0` or `1` based  intervals.
    """
    def makewindows(_c_len, index_start, windowsize):
        """create a file of the specified windows"""
        start = index_start
        end = min(_c_len, windowsize)
        starts = [start]
        ends = [end]
        while end < _c_len:
            end = min(end + windowsize, _c_len)
            ends.append(end)
            start += windowsize
            starts.append(start)
        return starts, ends

    if input.lower().endswith("fai"):
        with open(input, "r", encoding="utf-8") as fai, open(output, "w") as fout:
            for line in fai:
                lsplit = line.split("\t")
                contig = lsplit[0]
                c_len = int(lsplit[1])
                c_len += mode
                starts,ends = makewindows(c_len, mode, window)
                for startpos,endpos in zip(starts, ends):
                    fout.write(f"{contig}\t{startpos}\t{endpos}\n")
        return

    with pysam.FastxFile(input, persist = False) as FA, open(output, "w") as fout:
        for record in FA:
            chrom_name = record.name
            chom_len = len(record.sequence)
            starts,ends = makewindows(chom_len, mode, window)
            for startpos,endpos in zip(starts, ends):
                fout.write(f"{chrom_name}\t{startpos}\t{endpos}\n")
