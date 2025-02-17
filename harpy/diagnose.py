import os
import sys
import yaml
import subprocess
import rich_click as click
from rich.console import Console
from rich import print as rprint
from ._printing import print_error

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False))
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def diagnose(directory):
    """
    Run the Snakemake debugger to identify hang-ups 
    """
    directory = directory.rstrip("/")
    if not os.path.exists(f"{directory}/workflow/config.yaml"):
        print_error("missing config file", f"Target directory [blue]{directory}[/blue] does not contain the file [bold]workflow/config.yaml[/bold]")
        sys.exit(1)
    with open(f"{directory}/workflow/config.yaml", 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
    
    workflow = harpy_config["workflow"].replace(" ", "_")
    command = harpy_config["workflow_call"] + " --dry-run --debug-dag"
    console = Console(stderr=True)
    console.rule("[bold]Diagnosing Snakemake Job Graph", style = "green")
    try:
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        while True:
            output = process.stdout.readline()  # Read a line from stdout
            if not output and process.poll():  # Check if process is done
                break
            if output:  # Print the output if it's not empty
                if output.startswith("This was a dry-run"):
                    process.terminate()
                    exit(0)
                else:
                    rprint(f"[yellow]{output}", end = "", file= sys.stderr)
    except KeyboardInterrupt:
        console.print("")
        console.rule("[bold]Terminating Snakemake", style = "yellow")
        process.terminate()
        process.wait()
        sys.exit(1)