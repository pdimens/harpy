"""Harpy module to create HPC profiles for running snakemake on a cluster"""

import os
import rich_click as click
from .printfunctions import print_notice

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
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
            "commands": ["slurm", "htcondor", "lsf", "googlelifesci", "generic"],
        },
    ]
}

@click.command()
def generic():
    """Configuration for a generic scheduler"""
    outfile = "hpc/generic/config.yaml"
    os.makedirs("hpc/generic", exist_ok=True)
    if os.path.exists(outfile):
        print_notice(f"[blue bold]{outfile}[/blue bold] exists, overwriting.")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: cluster-generic\n")
        yml.write("default-resources:\n")
        yml.write("\tmem_mb: attempt * 2000\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 60\n")
        yml.write("retries: 1\n")
        yml.write("# command for submitting jobs\n")
        yml.write("cluster-generic-submit-cmd VALUE\n")
        yml.write("# [optional] command for retrieving job status\n")
        yml.write("cluster-generic-status-cmd VALUE\n")
        yml.write("# [optional] command for cancelling jobs-- expected to take one or more jobids as arguments\n")
        yml.write("cluster-generic-cancel-cmd VALUE\n")
        yml.write("# [optional] number of jobids to pass to cancel_cmd and if more are given, cancel_cmd will be called multiple times\n")
        yml.write("cluster-generic-cancel-nargs 20\n")
        yml.write("# [optional] command for sidecar process\n")
        yml.write("cluster-generic-sidecar-cmd VALUE\n")
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
    outfile = "hpc/lsf/config.yaml"
    os.makedirs("hpc/lsf", exist_ok=True)
    if os.path.exists(outfile):
        print_notice(f"[blue bold]{outfile}[/blue bold] exists, overwriting.")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: lsf\n")
        yml.write("default-resources:\n")
        yml.write("\tlsf_queue: \n")
        yml.write("\twalltime: 60 # minutes per job\n")
        yml.write("\tmem_mb: attempt * 2000\n")
        yml.write("\t# other args to pass to bsub\n")
        yml.write("\tlsf_extra: \n")
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
        yml.write("executor: htcondor\n")
        yml.write("default-resources:\n")
        yml.write("\tgetenv: True\n")
        yml.write("\trank: \n")
        yml.write("\tmax-retries: 1\n")
        yml.write("\trequest_memory: attempt * 2000\n")
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
def slurm():
    """Configuration for SLURM"""
    os.makedirs("hpc/slurm", exist_ok=True)
    outfile = "hpc/slurm/config.yaml"
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
hpc.add_command(generic)
hpc.add_command(googlelifesci)
