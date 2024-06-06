"""Harpy module to create HPC profiles for running snakemake on a cluster"""

import os
import rich_click as click
from .printfunctions import print_notice

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
def hpc():
    """
    Profile templates for cluster job submissions

    If running Harpy on an HPC system (or cluster), you can leverage Snakemake
    to handle all the job submissions on your behalf. This command creates templates
    for common HPC schedulers so you can run Harpy on a cluster with minimal friction.
    The subcommands create a `config.yml` in an `hpc/` directory.
    """

docstring = {
    "harpy hpc": [
        {
            "name": "Submission Systems",
            "commands": ["slurm", "htcondor", "lsf", "googlelifesci"],
        },
    ]
}

@click.command()
def slurm():
    """Configuration for SLURM"""
    outfile = "hpc/slurm/config.yaml"
    os.makedirs("hpc/slurm", exist_ok=True)
    if os.path.exists(outfile):
        print_notice(f"[blue bold]{outfile}[/blue bold] exists, overwriting.")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: slurm\n")
        yml.write("default-resources:\n")
        yml.write("\tslurm_account: $USER\n")
        yml.write("\tslurm_partition: regular\n")
        yml.write("\tmem_mb: attempt * 2000\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 60\n")
        yml.write("retries: 1\n")
        yml.write("\n# This section is for advanced copying into a scratch directory #\n")
        yml.write("#default-storage-provider: fs\n")
        yml.write("#local-storage-prefix: /home2/$USER\n")
        yml.write("#shared-fs-usage:\n")
        yml.write("#- persistence\n")
        yml.write("#- software-deployment\n")
        yml.write("#- sources\n")
        yml.write("#- source-cache\n")

@click.command()
def htcondor():
    """Configuration for HTCondor"""
    os.makedirs("hpc/htcondor", exist_ok=True)
    outfile = "hpc/htcondor/config.yaml"
    if os.path.exists(outfile):
        print_notice(f"[blue bold]{outfile}[/blue bold] exists, overwriting.")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: slurm\n")
        yml.write("default-resources:\n")
        yml.write("\tslurm_account: $USER\n")
        yml.write("\tslurm_partition: regular\n")
        yml.write("\tmem_mb: attempt * 2000\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 60\n")
        yml.write("retries: 1\n")
        yml.write("\n# This section is for advanced copying into a scratch directory #\n")
        yml.write("#default-storage-provider: fs\n")
        yml.write("#local-storage-prefix: /home2/$USER\n")
        yml.write("#shared-fs-usage:\n")
        yml.write("#- persistence\n")
        yml.write("#- software-deployment\n")
        yml.write("#- sources\n")
        yml.write("#- source-cache\n")

@click.command()
def lsf():
    """Configuration for LSF"""
    os.makedirs("hpc/lsf", exist_ok=True)
    outfile = "hpc/lsf/config.yaml"
    if os.path.exists(outfile):
        print_notice(f"[blue bold]{outfile}[/blue bold] exists, overwriting.")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: slurm\n")
        yml.write("default-resources:\n")
        yml.write("\tslurm_account: $USER\n")
        yml.write("\tslurm_partition: regular\n")
        yml.write("\tmem_mb: attempt * 2000\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 60\n")
        yml.write("retries: 1\n")
        yml.write("\n# This section is for advanced copying into a scratch directory #\n")
        yml.write("#default-storage-provider: fs\n")
        yml.write("#local-storage-prefix: /home2/$USER\n")
        yml.write("#shared-fs-usage:\n")
        yml.write("#- persistence\n")
        yml.write("#- software-deployment\n")
        yml.write("#- sources\n")
        yml.write("#- source-cache\n")

@click.command()
def googlelifesci():
    """Configuration for Google Life Sciences"""
    os.makedirs("hpc/googlelifesci", exist_ok=True)
    outfile = "hpc/googlelifesci/config.yaml"
    if os.path.exists(outfile):
        print_notice(f"[blue bold]{outfile}[/blue bold] exists, overwriting.")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("google-lifesciences: True\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 45\n")
        yml.write("retries: 1\n")
        yml.write("default-resources:\n")
        yml.write("\tmem_mb: attempt * 4000\n")
        yml.write("\tmem_mb_reduced: (attempt * 2000) * 0.9\n")
        yml.write("\tmachine_type: \"n2-standard-4\"\n")

hpc.add_command(slurm)
hpc.add_command(htcondor)
hpc.add_command(lsf)
hpc.add_command(googlelifesci)
