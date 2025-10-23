"""Harpy module to create HPC profiles for running snakemake on a cluster"""

import rich_click as click
from harpy.common.file_ops import fetch_template
from harpy.common.system_ops import package_absent

@click.command(panel = "HPC Configurations")
def hpc_generic():
    """
    Create a template config for a generic scheduler
    
    This command creates a configuration for a generic HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-cluster-generic` for the HPC job submission to work.
    """
    fetch_template("hpc-generic.yaml")
    package_absent("snakemake-executor-plugin-cluster-generic")

#def hpc_generic2():
#    pass
    #executor: cluster-generic
    #cluster-generic-submit-cmd:
    #  mkdir -p results/slurm_logs/{rule} &&
    #  sbatch
    #    --partition=regular
    #    --account=username
    #    --cpus-per-task={threads}
    #    --mem-per-cpu={resources.mem_per_cpu}
    #    --time={resources.time}
    #    --job-name=harpy-{rule}-{wildcards}
    #    --output=results/slurm_logs/{rule}/{rule}-%j-{wildcards}.out
    #    --error=results/slurm_logs/{rule}/{rule}-%j-{wildcards}.err
    #    --parsable
    #cluster-generic-status-cmd: status-sacct-robust.sh
    #cluster-generic-cancel-cmd: scancel
    #cluster-generic-cancel-nargs: 400
    #default-resources:
    #  - time="12:00:00"
    #  - mem_per_cpu=3200
    #  - tmpdir="/tmp"
    #restart-times: 2
    #max-jobs-per-second: 10
    #max-status-checks-per-second: 2
    #local-cores: 1
    #latency-wait: 60
    #cores: 800
    #jobs: 500
    #keep-going: True
    #rerun-incomplete: True

@click.command(panel = "HPC Configurations")
def hpc_lsf():
    """
    Create a template config for LSF
    
    This command creates a configuration for the LSF HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-lsf` for the HPC job submission to work.
    """
    fetch_template("hpc-lsf.yaml")

    package_absent("snakemake-executor-plugin-lsf")

@click.command(panel = "HPC Configurations")
def hpc_slurm():
    """
    Create a template config for SLURM
    
    This command creates a configuration for the SLURM HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-slurm` for the HPC job submission to work.
    """
    fetch_template("hpc-slurm.yaml")
    package_absent("snakemake-executor-plugin-slurm")

@click.command(panel = "HPC Configurations")
def hpc_googlebatch():
    """
    Create a template config for Google Batch
    
    This command creates a configuration for the Google Batch scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-googlebatch` for the HPC job submission to work.
    """
    fetch_template("hpc-googlebatch.yaml")
    package_absent("snakemake-executor-plugin-googlebatch")
