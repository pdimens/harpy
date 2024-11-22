"""View the latest log or snakefile of a workflow"""

import os
import sys
import glob
import subprocess
import rich_click as click
from ._printing import print_error
from ._validations import is_gzip

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False))
@click.option('-s', '--snakefile',  is_flag = True, show_default = True, default = False, help = "View the snakefile, not the log file")
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def view(directory, snakefile):
    """
    View a workflow log file or snakefile

    This convenience command lets you view the latest workflow log file
    of a Harpy output directory. You can use `--snakefile` to view the workflow
    snakefile instead. The output will be printed to your screen via the `less` command,
    so use standard `less` keyboard shortcuts to navigate the output, e.g.:
    
    | key                     | function                   |
    | :---------------------- | :------------------------- |
    | `up arrow`/`down arrow` | scroll up/down             |
    | `Page Up`/`Page Down`   | scroll up/down, but faster |
    | `q`                     | exit                       |
    """
    # check if there is a workflow or log folder
    # and whether the expected files are in there
    err = 0
    if snakefile:
        files = [i for i in glob.iglob(f"{directory}/workflow/*.smk")]
        err_dir = f"{directory}/workflow/"
        err_file = "snakefiles"
        if not os.path.exists(f"{directory}/workflow"):
            err = 1
        elif not files:
            err = 2
    else:
        files = [i for i in glob.iglob(f"{directory}/logs/snakemake/*.log*")]
        err_dir = f"{directory}/logs/snakemake/"
        err_file = "log files"
        if not os.path.exists(f"{directory}/logs/snakemake"):
            err = 1
        elif not files:
            err = 2
    if err == 1:
        print_error(
            "Directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/blue], but that directory was not found. Please check that this is the correct folder."
        )
        sys.exit(1)
    elif err == 2:
        print_error(
            "File not found", 
            f"There are no {err_file} in [blue]{err_dir}[/blue]. Please check that this is the correct folder."
        )
        sys.exit(1)
    # sort and pull only the most recent file (based on modification time)
    file = sorted(files, key = os.path.getmtime)[-1]
    if not os.access(file, os.R_OK):
        print_error(
            "Incorrect permissions",
            f"[blue]{file}[/blue] does not have read access. Please check the file permissions."
        )
        sys.exit(1)

    cat_cmd = "zcat" if is_gzip(file) else "cat"
    stream = subprocess.Popen([cat_cmd, file], stdout=subprocess.PIPE)
    pygment = subprocess.Popen(["pygmentize", "-l", "yaml"], stdin = stream.stdout, stdout = subprocess.PIPE)
    subprocess.run(["less", "-R"], stdin = pygment.stdout) 
