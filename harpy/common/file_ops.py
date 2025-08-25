"""Module with helper function to set up Harpy workflows"""

import os
import glob
import gzip
import shutil
from pathlib import Path
import importlib.resources as resources
from harpy.common.printing import print_error

def fetch_snakefile(workdir: str, target: str) -> None:
    """
    Retrieve the target harpy rule and write it into the workdir as workflow.smk
    """
    os.makedirs(workdir, exist_ok= True)
    dest_file = os.path.join(workdir,"workflow.smk")
    source_file = resources.files("harpy.snakefiles") / target
    try:
        with resources.as_file(source_file) as _source:
            shutil.copy2(_source, dest_file)
    except (FileNotFoundError, KeyError):
        print_error(
            "snakefile missing",
            f"The required snakefile [blue bold]{target}[/] was not found in the Harpy installation.",
            "There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be an issue with your conda/mamba environment or configuration."
        )

def filepath(infile: str) -> str:
    """returns a posix-formatted absolute path of infile"""
    return Path(infile).resolve().as_posix()

def symlink(original: str, destination: str) -> None:
    """Create a symbolic link from original -> destination if the destination doesn't already exist."""
    if not (Path(destination).is_symlink() or Path(destination).exists()):
        Path(destination).symlink_to(Path(original).resolve())

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
