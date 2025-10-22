import glob
import os
import sys
import yaml
import subprocess
import rich_click as click
from harpy.common.printing import print_error, CONSOLE, print_shellcmd_simple
from harpy.common.file_ops import safe_read

@click.group(options_metavar='', context_settings={"help_option_names" : []})
def diagnose():
    """
    Attempt to resolve workflow errors
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def stall(directory):
    """
    Run the Snakemake debugger to identify why a workflow stalled

    This will run Snakemake with the `--dry-run` and `--debug-dag` options,
    printing the diagnostics to the terminal.
    """
    directory = directory.rstrip("/")
    PROFILE_FILE = os.path.join(directory, "workflow", "config.yaml")
    CONFIG_FILE = os.path.join(directory, "workflow", "workflow.yaml")

    if not os.path.exists(CONFIG_FILE):
        print_error("missing workflow config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/workflow.yaml[/]")
    if not os.path.exists(PROFILE_FILE):
        print_error("missing snakemake config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/config.yaml[/]")

    with open(CONFIG_FILE, 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
    
    command = harpy_config["snakemake"]["absolute"]
    # prefix the new arguments, in case a positional argument was added at the end by user
    command = command.replace("snakemake -", "snakemake --sdm env-modules --dry-run --debug-dag -")

    CONSOLE.rule("[bold]Diagnosing Snakemake Job Graph", style = "green")
    try:
        process = subprocess.Popen(
            command.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8'
        )
        while True:
            output = process.stdout.readline()
            error = process.stderr.readline()
            if not output and not error and process.poll() is not None:
                break
            if error:
                CONSOLE.print(error, end="", style= "red")
                # error usually prints more than one line, so this will make sure all
                # consecutive stderr text will be printed together
                while error:
                    error = process.stderr.readline()
                    CONSOLE.print(error, end="", style = "red")
            if output:
                if output.startswith("This was a dry-run"):
                    process.terminate()
                    exit(0)
                if "Exception" in output:
                    while output:
                        CONSOLE.print(output, end = "")
                        output = process.stdout.readline()
                elif output.lstrip().startswith("["):
                    CONSOLE.print(f"\n{output}", end = "", highlight=False, style = "blue")
                else:
                    CONSOLE.print(output, end="", style = "yellow")
    except Exception as e:
        CONSOLE.print("")
        CONSOLE.rule("[bold]End of diagnosis", style = "yellow")
        process.terminate()
        process.wait()
        sys.exit(1)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def snakemake(directory):
    """
    Run the Snakemake workflow directly without Harpy intervention

    This will run Snakemake without any of the convenience features,
    which can sometimes reveal errors that weren't captured by Harpy. Requries
    the `workflow/config.yaml` and `workflow/workflow.yaml` files to be in
    `DIRECTORY`. 
    """
    directory = directory.rstrip("/")
    PROFILE_FILE = os.path.join(directory, "workflow", "config.yaml")
    CONFIG_FILE = os.path.join(directory, "workflow", "workflow.yaml")

    if not os.path.exists(CONFIG_FILE):
        print_error("missing workflow config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/workflow.yaml[/]")
    if not os.path.exists(PROFILE_FILE):
        print_error("missing snakemake config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/config.yaml[/]")

    with open(CONFIG_FILE, 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
    
    command = harpy_config["snakemake"]["absolute"]

    os.system(command)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def rule(directory):
    """
    Directly run the first rule that caused the workflow failure

    The rule is identified in the most recent Snakemake log file in
    `DIRECTORY` as the first one with the text `Error in rule ____`.
    """
    directory = directory.rstrip("/")
    if not os.path.exists(f'{directory}/logs/snakemake/'):
        print_error("missing log folder", f"Target directory [blue]{directory}[/] does not contain the folder [bold]logs/snakemake[/]")
    # get the lastest snakemake log file
    list_of_files = glob.glob(f'{directory}/logs/snakemake/*')
    latest_log = max(list_of_files, key=os.path.getctime)

    CONSOLE.rule(f"Latest log: {os.path.basename(latest_log)}", style = "yellow")
    _found = False
    _shellblock = False
    conda = ""
    container = ""
    cmd = []
    with safe_read(latest_log) as logfile:
        for line in logfile:
            if "Error in rule" in line:
                CONSOLE.print(f"Failing rule: {line.strip().removeprefix('Error in rule').rstrip(':')}", style = "yellow")
                _found = True
                continue
            if _found:
                if "conda-env:" in line:
                    conda = "source activate " + line.strip().split()[-1]
                    continue
                if "container:" in line:
                    container = line.strip().split()[-1]
                if "shell:" in line:
                    _shellblock = True
                    _ = line.strip().replace("shell:","").split()
                    if _:
                        cmd.append(_)
                    continue
            if _shellblock:
                if line.strip() == "(command exited with non-zero exit code)":
                    break
                cmd.append(line.strip())

    if conda:
        print_shellcmd_simple("\n".join(cmd))
        os.system("\n".join([conda, f"cd {directory}", *cmd]))
    elif container:
        print_shellcmd_simple(f"""
apptainer exec {container} bash -c '
{"\n".join([*cmd])}
'
"""
        )
        os.system(f"""
cd {directory}
apptainer exec {container} bash -c '
{"\n".join([*cmd])}
'
"""
        )
    else:
        print_shellcmd_simple("\n".join(cmd))
        os.system("\n".join([f"cd {directory}", *cmd]))

diagnose.add_command(stall)
diagnose.add_command(snakemake)
diagnose.add_command(rule)