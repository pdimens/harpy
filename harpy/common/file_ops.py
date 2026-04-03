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
import curses
from datetime import datetime
from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import get_formatter_by_name
from click import echo_via_pager
from rich.console import Console
from rich.prompt import Prompt
from rich.syntax import Syntax
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

#========== harpy view ==============#

def check_terminal_colors():
    # Initialize curses and always tear down
    try:
        _ = curses.initscr()
        # Check if the terminal supports colors
        if not curses.has_colors():
            return 0
        curses.start_color()
        num_colors = curses.COLORS
        # Determine the color type based on the number of colors
        if num_colors <= 8:
            return 8
        else:
            return 256
    except curses.error:
        # Non-interactive/unsupported terminals
        return 0
    finally:
        try:
            curses.endwin()
        except curses.error:
            pass

def parse_file(infile: str):
    '''
    Print file contents via pygmentized `less`.
    '''
    hp = HarpyPrint()
    if not os.access(infile, os.R_OK):
        hp.error(
            "incorrect permissions",
            f"[blue]{infile}[/] does not have read access. Please check the file permissions."
        )
    n_colors = check_terminal_colors()
    if n_colors <= 8:
        formatter = get_formatter_by_name("terminal")
    else:
        formatter = get_formatter_by_name("terminal256")

    def _read_file(x: str):
        compressed = is_gzip(x)
        opener = gzip.open if compressed else open
        mode = "rt" if compressed else "r"
        lexer = get_lexer_by_name("yaml")
        with opener(x, mode) as f:
            for line in f:
                yield highlight(line, lexer, formatter)
    os.environ["PAGER"] = "less -R"
    echo_via_pager(_read_file(infile), color = n_colors > 0)

def parse_error(infile: str):
    '''
    Print syntax-highlighted error in snakemake log.
    '''
    hp = HarpyPrint()
    if not os.access(infile, os.R_OK):
        hp.error(
            "incorrect permissions",
            f"[blue]{infile}[/] does not have read access. Please check the file permissions."
        )
    with safe_read(infile) as f:
        result = ""
        for line in f:
            if "error" in line.lower() or "exception" in line.lower():
                result += line
                while True:
                    _line = f.readline()
                    if not _line:
                        break
                    result += _line
                break
        hp.rule(f"[default]{os.path.relpath(infile)}", style = 'blue')
        hp.print(Syntax(result, lexer= "yaml", background_color= 'default'), soft_wrap=True)

def choose_logfile(directory: str, choose:bool) -> str:
    hp = HarpyPrint()
    err_dir = os.path.join(directory, ".snakemake", "log")
    err_file = "There are no log files"
    if not os.path.exists(err_dir):
        hp.error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that this is the correct folder."
        )
    files = [i for i in glob.iglob(f"{err_dir}/*.log*")]        
    if not files:
        hp.error(
            "files not found", 
            f"{err_file} in [blue]{err_dir}[/]. Please check that this is the correct folder."
        )

    files = sorted(files, key = os.path.getmtime, reverse = True)
    if choose and len(files) > 1:
        console = Console()
        console.print()
        #console.rule('Snakemake Log Files', style = "green")
        _tb = hp.table()
        _tb.show_header=True
        _tb.add_column("[bold green]#", style="bold green", min_width=2)
        _tb.add_column("[dim yellow]Last Modification",style = "dim yellow")
        _tb.add_column("Log File", justify="right", no_wrap=True)
        for i,j in enumerate(files,1):
            filename = os.path.basename(j).removesuffix(".snakemake.log") + "[dim].snakemake.log[/]"
            modtime = datetime.fromtimestamp(os.path.getmtime(j)).strftime('%Y-%m-%d %H:%M')
            _tb.add_row(str(i), modtime , filename)
        console.print(_tb)
        selection = Prompt.ask(
            "\n[bold blue]Select a log file by number ([bold green]#[/])[/]",
            choices=list(str(i) for i in range(1,len(files) + 1)),
            show_choices=False
        )
        
        selected_idx = int(selection) - 1
        target_file = files[selected_idx]
    else:
        target_file = files[0]
    return target_file