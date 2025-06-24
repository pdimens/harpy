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
    if not os.path.exists(f"{directory}/workflow/workflow.yaml"):
        print_error("missing workflow config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/workflow.yaml[/]")
        sys.exit(1)
    if not os.path.exists(f"{directory}/workflow/config.yaml"):
        print_error("missing snakemake config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/config.yaml[/]")
        sys.exit(1)

    with open(f"{directory}/workflow/workflow.yaml", 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
    
    command = harpy_config["snakemake"]["absolute"]
    # prefix the new arguments, in case a positional argument was added at the end by user
    command = command.replace("snakemake -", "snakemake --dry-run --debug-dag -")
    console = Console(stderr=True)
    console.rule("[bold]Diagnosing Snakemake Job Graph", style = "green")
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
                rprint(f"[red]{error}", end="", file=sys.stderr)
                # error usually prints more than one line, so this will make sure all
                # consecutive stderr text will be printed together
                while error:
                    error = process.stderr.readline()
                    rprint(f"[red]{error}", end="", file=sys.stderr)
            if output:
                if output.startswith("This was a dry-run"):
                    process.terminate()
                    exit(0)
                else:
                    rprint(f"[yellow]{output}", end="", file=sys.stderr)
    except:
        console.print("")
        console.rule("[bold]End of diagnosis", style = "yellow")
        process.terminate()
        process.wait()
        sys.exit(1)