import glob
import os
import sys
import yaml
import subprocess
import rich_click as click
from harpy.common.printing import HarpyPrint
from harpy.common.file_ops import safe_read

@click.group(options_metavar='')
@click.help_option('--help', hidden = True)
def diagnose():
    """
    Attempt to resolve workflow errors
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.help_option('--help', hidden = True)
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def stall(directory):
    """
    Run the Snakemake debugger to identify why a workflow stalled

    This will run Snakemake with the `--dry-run` and `--debug-dag` options,
    printing the diagnostics to the terminal.
    """
    hp = HarpyPrint()
    directory = directory.rstrip("/")
    PROFILE_FILE = os.path.join(directory, "workflow", "profile.yaml")
    CONFIG_FILE = os.path.join(directory, "workflow", "workflow.yaml")

    if not os.path.exists(CONFIG_FILE):
        hp.error("missing workflow config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/workflow.yaml[/]")
    if not os.path.exists(PROFILE_FILE):
        hp.error("missing snakemake config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/profile.yaml[/]")

    with open(CONFIG_FILE, 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
    
    command = harpy_config["snakemake"]["absolute"]
    # prefix the new arguments, in case a positional argument was added at the end by user
    command = command.replace("snakemake -", "snakemake --sdm env-modules --dry-run --debug-dag -")

    hp.rule("[bold]Diagnosing Snakemake Job Graph", style = "green")
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
                hp.print(error, end="", style= "red")
                # error usually prints more than one line, so this will make sure all
                # consecutive stderr text will be printed together
                while error:
                    error = process.stderr.readline()
                    hp.print(error, end="", style = "red")
            if output:
                if output.startswith("This was a dry-run"):
                    process.terminate()
                    exit(0)
                if "Exception" in output:
                    while output:
                        hp.print(output, end = "")
                        output = process.stdout.readline()
                elif output.lstrip().startswith("["):
                    hp.print(f"\n{output}", end = "", highlight=False, style = "blue")
                else:
                    hp.print(output, end="", style = "yellow")
    except Exception as e:
        hp.print("")
        hp.rule("[bold]End of diagnosis", style = "yellow")
        process.terminate()
        process.wait()
        sys.exit(1)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.help_option('--help', hidden = True)
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def rule(directory):
    """
    Directly run the first rule that caused the workflow failure

    The rule is identified in the most recent Snakemake log file in
    `DIRECTORY/logs/snakemake` as the first one with the text `Error in rule ____`. 
    If the failing rule is missing inputs (e.g. they were temporary), Harpy will run Snakemake first to generate
    those files, then execute the failing rule directly (i.e. without Snakemake).
    """
    hp = HarpyPrint()
    directory = directory.rstrip("/")
    PROFILE_FILE = os.path.join(directory, "workflow", "profile.yaml")
    CONFIG_FILE = os.path.join(directory, "workflow", "workflow.yaml")

    if not os.path.exists(f'{directory}/.snakemake/log'):
        hp.error("missing log folder", f"Target directory [blue]{directory}[/] does not contain the folder [bold].snakemake/log[/]")
    # get the lastest snakemake log file
    list_of_files = glob.glob(f'{directory}/.snakemake/log/*')
    if not list_of_files:
        hp.error("missing log files", f"Log directory [blue]{directory}/.snakemake/log[/] does not have any log files in it")

    latest_log = max(list_of_files, key=os.path.getctime)

    hp.rule(f"Latest log: [bold default]{os.path.basename(latest_log)}", style = "yellow")
    failed_rule = ""
    _shellblock = False
    infiles = []
    conda = ""
    container = ""
    cmd = []
    with safe_read(latest_log) as logfile:
        for line in logfile:
            if "Error in rule" in line:
                failed_rule += line.strip().removeprefix('Error in rule').rstrip(':')
                continue
            if failed_rule:
                if "input:" in line:
                    for i in line.strip().replace("input: ", "").split(", "):
                        if not os.path.exists(os.path.join(directory,i)):
                            infiles.append(i)
                if "conda-env:" in line:
                    conda = "source activate " + line.strip().split()[-1]
                    continue
                if "container:" in line:
                    container = line.strip().split()[-1]
                if "shell:" in line:
                    _shellblock = True
                    _ = line.strip().replace("shell:","").strip().split()
                    if _:
                        cmd.append(_)
                    continue
            if _shellblock:
                if line.strip() == "(command exited with non-zero exit code)":
                    break
                cmd.append(line.strip())
    if failed_rule:
        hp.log(f"Failing rule: [yellow]{failed_rule.lstrip()}", newline=True)
    else:
        hp.log(f"No errors found in {os.path.basename(latest_log)}", style = "green", markup=False, highlight=False)
        sys.exit(0)
    if infiles:
        if not os.path.exists(CONFIG_FILE):
            hp.error("missing workflow config", f"The failing rule is missing inputs, which requires Snakemake to be re-run so they can be generated, but target directory [blue]{directory}[/] does not contain the file [bold]workflow/workflow.yaml[/]")
        if not os.path.exists(PROFILE_FILE):
            hp.error("missing snakemake config", f"The failing rule is missing inputs, which requires Snakemake to be re-run so they can be generated, but target directory [blue]{directory}[/] does not contain the file [bold]workflow/profile.yaml[/]")
        hp.log("Missing input files:\n  [yellow]" + '\n  '.join(infiles))

        with open(CONFIG_FILE, 'r', encoding="utf-8") as f:
            harpy_config = yaml.full_load(f)
            command = harpy_config["snakemake"]["absolute"]

        command += f" --quiet --no-temp {' '.join(infiles)}"
        hp.log("Rerunning Snakemake to generate inputs")
        hp.shell(command)
        sm = os.system(command)
        if sm != 0:
            hp.error("workflow error", "Harpy attempted to regenerate the input files necessary to run the failed rule directly, but that seemed to fail too. You may want to try manually rerunning the step(s) that failed.")
    else:
        hp.log("Missing input files: [green]None")
    hp.log("Running failed code block")
    if conda:
        #hp.shell("\n".join(cmd))
        hp.shell("\n".join([conda, f"cd {directory}", *cmd]), rules = True)

        os.system("\n".join([conda, f"cd {directory}", *cmd]))
    elif container:
        joined_cmd = "\n".join([*cmd])
        hp.shell(f"""apptainer exec {container} bash -c '\n{joined_cmd}\n'""", rules = True)
        os.system(f"""
cd {directory}
apptainer exec {container} bash -c '
{joined_cmd}
'
"""
        )
    else:
        hp.shell("\n".join([f"cd {directory}", *cmd]), rules=True)
        os.system("\n".join([f"cd {directory}", *cmd]))

diagnose.add_command(stall)
diagnose.add_command(rule)