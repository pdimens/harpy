from .helperfunctions import get_samples_from_fastq, print_error, print_solution

import os
import sys
import rich_click as click

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-o', '--output', required = True, type=str, metavar = "Output file", help = 'Name of output file')
@click.option('-s', '--system', required = True, type = click.Choice(["slurm"], case_sensitive = False), help = 'HPC scheduling system [slurm]')
def hpc(output, system):
    """
    Create a config file to run Harpy on an HPC

    With this command you can generate a HPC profile for running
    Harpy on a cluster.
    """
    if hpc == "sge":
        click.echo("Unfortunately, SGE is not yet supported.")
        exit(1)
    if hpc == 'slurm':
        if os.path.exists(output):
            overwrite = input(f"File {output} already exists, overwrite (no|yes)?  ").lower()
            if overwrite not in ["yes", "y"]:
                click.echo("Please suggest a different name for the output file")
                exit(0)
        with open(output, "w") as yml:
            yml.write(
                """cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --partition={resources.partition}
    --qos={resources.qos}
    --cpus-per-task={threads}
    --mem={resources.mem_mb}
    --job-name=smk-{rule}-{wildcards}
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
default-resources:
  - partition=<name-of-default-partition>
  - qos=<name-of-quality-of-service>
  - mem_mb=1000
restart-times: 3
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 500
keep-going: True
rerun-incomplete: True
printshellcmds: False
scheduler: greedy
use-conda: False
"""
        )
    click.echo(f"Created HPC profile {output}. Replace the \'partition\' and \'qos\' placeholders in {output} with options relevant to your system.", file = sys.stderr)