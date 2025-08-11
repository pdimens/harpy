import os
import sys
import yaml
import subprocess
import rich_click as click
from harpy.common.printing import print_error, CONSOLE

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False))
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def diagnose(directory):
    """
    Run the Snakemake debugger to identify hang-ups 
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
    command = command.replace("snakemake -", "snakemake --dry-run --debug-dag -")
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
                else:
                    CONSOLE.print(output, end="", style = "yellow")
    except:
        CONSOLE.print("")
        CONSOLE.rule("[bold]End of diagnosis", style = "yellow")
        process.terminate()
        process.wait()
        sys.exit(1)