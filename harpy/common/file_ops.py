"""Module with helper function to set up Harpy workflows"""

import os
import glob
import gzip
import importlib.resources as resources
import shutil
import sys
from pathlib import Path
from harpy.common.printing import print_error

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
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        _out = open(outfile, "w")
    else:
        _out = sys.stdout
    try:
        with resources.as_file(source_file) as _source, open(_source, 'r') as f:
            _out.write(f.read() + "\n")
    except (FileNotFoundError, KeyError):
        print_error(
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
