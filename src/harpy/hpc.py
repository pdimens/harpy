from .printfunctions import print_error, print_solution
from .fileparsers import get_samples_from_fastq
import os
import sys
import rich_click as click

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/ ")
@click.option('-o', '--output', required = True, type=str, metavar = "Output file", help = 'Name of output file')
@click.option('-s', '--system', required = True, type = click.Choice(["slurm"], case_sensitive = False), help = 'HPC scheduling system [slurm]')
def hpc(output, system):
    """
    Create a config file to run Harpy on an HPC

    With this command you can generate a HPC profile for running
    Harpy on a cluster.
    """
    if os.path.exists(output):
        overwrite = input(f"File {output} already exists, overwrite (no|yes)?  ").lower()
        if overwrite not in ["yes", "y"]:
            click.echo("Please suggest a different name for the output file")
            exit(0)
    with open(output, "w") as outfile:
        outfile.write("cluster:\n")
        outfile.write("\tmkdir -p logs/{rule} &&\n")
        outfile.write("\tsbatch\n")
        outfile.write("\t\t--partition={resources.partition}\n")
        outfile.write("\t\t--qos={resources.qos}\n")
        outfile.write("\t\t--cpus-per-task={threads}\n")
        outfile.write("\t\t--mem={resources.mem_mb}\n")
        outfile.write("\t\t--job-name=smk-{rule}-{wildcards}\n")
        outfile.write("\t\t--output=logs/{rule}/{rule}-{wildcards}-%j.out\n")
        outfile.write("default-resources:\n")
        outfile.write("\t- partition=<name-of-default-partition>\n")
        outfile.write("\t- qos=<name-of-quality-of-service>\n")
        outfile.write("\t- mem_mb=1000\n")
        outfile.write("restart-times: 3\n")
        outfile.write("max-jobs-per-second: 10\n")
        outfile.write("max-status-checks-per-second: 1\n")
        outfile.write("local-cores: 1\n")
        outfile.write("latency-wait: 60\n")
        outfile.write("jobs: 500\n")
        outfile.write("keep-going: True\n")
        outfile.write("rerun-incomplete: True\n")
        outfile.write("printshellcmds: False\n")
        outfile.write("scheduler: greedy\n")
        outfile.write("use-conda: False\n")
    click.echo(f"Created HPC profile {output}. Replace the \'partition\' and \'qos\' placeholders in {output} with options relevant to your system.", file = sys.stderr)